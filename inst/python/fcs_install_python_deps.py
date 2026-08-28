#!/usr/bin/env python3
"""
Lightweight dependency checker/installer for FCSimple's Python helpers.

This script does NOT import heavy libraries (numba/umap/pynndescent) to avoid
triggering JIT compilation or native runtime initialization. It uses
``python -m pip show`` to detect installed distributions for the same Python
interpreter that executes the script. Optionally it can install missing
packages and optionally run the bundled precompile script in a separate
process.

Usage:
  python fcs_install_python_deps.py [--install] [--precompile] [--build-optsne] [--packages pkg1,pkg2]

Exit code: 0 = all packages present (or successfully installed), 1 = missing
packages (and not installed) or installation failure, 2 = precompile failure
"""
import argparse
import subprocess
import sys
import os

DEFAULT_PACKAGES = [
    "pandas",
    "numpy",
    "scipy",
    "scikit-learn",
    "umap-learn",
    "numba",
    "pynndescent",
    "tqdm",
    "cffi",
    "openTSNE",
    "igraph",
    "leidenalg",
    "flowkit"
]


def pip_show(python_exe, pkg):
    try:
        res = subprocess.run([python_exe, "-m", "pip", "show", pkg], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return res.returncode == 0
    except Exception:
        return False


def pip_install(python_exe, pkg):
    try:
        subprocess.check_call([python_exe, "-m", "pip", "install", pkg])
        return True
    except subprocess.CalledProcessError:
        return False


def _find_in_scripts(name):
    """Locate an executable in the Python Scripts directory (pip wheel output)."""
    import shutil
    import sysconfig

    found = shutil.which(name)
    if found:
        return found

    scripts = sysconfig.get_path("scripts")
    if scripts:
        candidate = os.path.join(scripts, name)
        if os.path.exists(candidate):
            return candidate

    return None


def _find_cmake():
    """Locate the cmake executable (PATH first, then pip wheel Scripts dir)."""
    return _find_in_scripts("cmake.exe" if os.name == "nt" else "cmake")


def _detect_cxx_compiler():
    """Return (compiler_path, compiler_bin_dir) for a usable C++ compiler.

    Windows supports MSYS2 (UCRT64/MINGW64/CLANG64) MinGW-w64, MSVC, or any
    compiler reachable via PATH. Returns (None, None) when nothing is found.
    """
    import shutil

    for name in ("g++", "gcc", "c++", "clang++"):
        p = shutil.which(name)
        if p:
            return p, os.path.dirname(p)

    if os.name == "nt":
        # MSYS2 under well-known roots; bin dirs are typically not on PATH.
        roots = []
        explicit = os.environ.get("MSYS2_ROOT")
        if explicit:
            roots.append(explicit)
        for base in (r"C:\msys64", r"D:\msys64"):
            if os.path.isdir(base) and base not in roots:
                roots.append(base)

        for root in roots:
            for env_sub in ("ucrt64", "mingw64", "clang64"):
                bin_dir = os.path.join(root, env_sub, "bin")
                for tool in ("g++.exe", "gcc.exe"):
                    candidate = os.path.join(bin_dir, tool)
                    if os.path.isfile(candidate):
                        return candidate, bin_dir

        # MSVC via the standard vswhere locator.
        vswhere = r"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"
        if os.path.isfile(vswhere):
            try:
                root = subprocess.check_output(
                    [vswhere, "-latest", "-products", "*", "-requires",
                     "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                     "-property", "installationPath"],
                    text=True,
                ).strip()
                vc_tools = os.path.join(root, "VC", "Tools", "MSVC")
                for dirpath, _dirs, files in os.walk(vc_tools):
                    if "cl.exe" in files:
                        cl = os.path.join(dirpath, "cl.exe")
                        return cl, dirpath
            except Exception:
                pass

    return None, None


_MINGW_RUNTIME_DLLS = (
    "libstdc++-6.dll",
    "libgcc_s_seh-1.dll",
    "libgcc_s_dw2-1.dll",
    "libwinpthread-1.dll",
    "libgomp-1.dll",
    "libatomic-1.dll",
    "libssp-0.dll",
    "libquadmath-0.dll",
)


def _copy_mingw_runtime_dlls(compiler_bin_dir, python_exe):
    """Copy MinGW runtime DLLs next to the built extension.

    A MinGW-built .dll depends on libstdc++/libgcc/libwinpthread/libgomp, whose
    DLLs live in the MSYS2 bin directory (not on the system PATH). Copying them
    into the installed MulticoreTSNE package directory lets the Python wrapper
    load the extension without PATH hacks.
    """
    import shutil
    import sysconfig

    try:
        pkg_dir = subprocess.check_output(
            [python_exe, "-c",
             "import MulticoreTSNE, os, sys; sys.stdout.write(os.path.dirname(os.path.abspath(MulticoreTSNE.__file__)))"],
            text=True,
        ).strip()
    except Exception:
        pkg_dir = os.path.join(sysconfig.get_paths().get("purelib", ""), "MulticoreTSNE")

    if not pkg_dir or not os.path.isdir(pkg_dir):
        return

    copied = []
    for dll in _MINGW_RUNTIME_DLLS:
        src = os.path.join(compiler_bin_dir, dll)
        dst = os.path.join(pkg_dir, dll)
        if os.path.isfile(src) and not os.path.exists(dst):
            try:
                shutil.copy2(src, dst)
                copied.append(dll)
            except OSError:
                pass

    if copied:
        print("Copied MinGW runtime DLLs next to extension:", ", ".join(copied))


def build_optsne(python_exe):
    """Build/install the bundled opt-SNE package (MulticoreTSNE).

    opt-SNE is not on PyPI and requires cmake + a C++ compiler. The source is
    vendored under inst/python/optsne so the build does not depend on the
    upstream GitHub repository. On Windows this locates the MSYS2/MinGW or
    MSVC compiler, installs cmake/ninja via pip when missing, and copies the
    MinGW runtime DLLs needed at import time.

    This is idempotent: if MulticoreTSNE is already installed, it succeeds
    without requiring cmake or a C++ compiler.
    """
    if pip_show(python_exe, "MulticoreTSNE"):
        print("opt-SNE (MulticoreTSNE) is already installed.")
        return True

    here = os.path.dirname(os.path.abspath(__file__))
    optsne_dir = os.path.join(here, "optsne")
    if not os.path.exists(os.path.join(optsne_dir, "setup.py")):
        print("Bundled opt-SNE source not found at", optsne_dir)
        return False

    compiler_path, compiler_dir = _detect_cxx_compiler()
    if compiler_path is None:
        print(
            "No C++ compiler found. On Windows install 'Microsoft C++ Build Tools' "
            "(Desktop development with C++ workload) or MSYS2 MinGW-w64."
        )
        return False
    print("Using C++ compiler:", compiler_path)

    # Assemble PATH additions and a CMake generator. Ninja is used when
    # available so the build does not depend on mingw32-make / make.
    path_additions = []
    generator = None

    cmake_exe = _find_cmake()
    if cmake_exe is None:
        print("cmake not found; installing it via pip ...")
        try:
            subprocess.check_call([python_exe, "-m", "pip", "install", "cmake"])
        except subprocess.CalledProcessError:
            print("Failed to install cmake via pip.")
        cmake_exe = _find_cmake()
    if cmake_exe is None:
        print(
            "cmake was not found on PATH. Install it (e.g. `pip install cmake` "
            "or the official installer) before building opt-SNE."
        )
        return False
    path_additions.append(os.path.dirname(cmake_exe))

    ninja_exe = _find_in_scripts("ninja.exe" if os.name == "nt" else "ninja")
    if ninja_exe is None and os.name == "nt":
        print("ninja not found; installing it via pip ...")
        try:
            subprocess.check_call([python_exe, "-m", "pip", "install", "ninja"])
        except subprocess.CalledProcessError:
            print("Failed to install ninja via pip.")
        ninja_exe = _find_in_scripts("ninja.exe")

    if ninja_exe is not None:
        path_additions.append(os.path.dirname(ninja_exe))
        generator = "Ninja"

    if compiler_dir:
        path_additions.append(compiler_dir)

    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(path_additions + [env.get("PATH", "")])
    if generator:
        env["CMAKE_GENERATOR"] = generator

    # Tell the CMake build recipe which compiler to use (essential for
    # MSYS2/MinGW/RTools toolchains that live outside the default PATH).
    # Match the tool name case-insensitively and regardless of a .exe suffix
    # (e.g. "g++.EXE" / "g++.exe" / "gcc").
    base = os.path.splitext(os.path.basename(compiler_path))[0].lower()
    if base in ("g++", "gcc", "clang++", "c++"):
        env["CMAKE_CXX_COMPILER"] = compiler_path

    print("Building opt-SNE from bundled source:", optsne_dir)
    try:
        subprocess.check_call([python_exe, "-m", "pip", "install", optsne_dir], env=env)
    except subprocess.CalledProcessError:
        print(
            "\nopt-SNE build failed. Ensure a C++ compiler (MSVC or MinGW-w64) and "
            "cmake are available, then retry."
        )
        return False

    # MinGW builds may need their runtime DLLs resolvable at import time (only
    # necessary for dynamically-linked MinGW/MSYS2 GCC; static RTools builds
    # have no such DLLs, so this becomes a no-op).
    if os.name == "nt" and compiler_path and "g++" in os.path.basename(compiler_path).lower():
        _copy_mingw_runtime_dlls(compiler_dir, python_exe)

    return True


def run_precompile(python_exe):
    here = os.path.dirname(os.path.abspath(__file__))
    precompile = os.path.join(here, "fcsimple_precompile_numba.py")
    if not os.path.exists(precompile):
        print("Precompile script not found; skipping precompile.")
        return True
    env = os.environ.copy()
    env["NUMBA_NUM_THREADS"] = "1"
    # run precompile in separate process so it can import numba/umap safely
    try:
        res = subprocess.run([python_exe, precompile], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print("Precompile failed:")
            print(res.stdout)
            print(res.stderr)
            return False
        return True
    except Exception as e:
        print("Precompile invocation error:", e)
        return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--install", action="store_true", help="Install missing packages via pip")
    p.add_argument("--precompile", action="store_true", help="Run one-time numba/umap precompile (single-threaded)")
    p.add_argument("--build-optsne", action="store_true", help="Build/install the bundled opt-SNE (MulticoreTSNE) package")
    p.add_argument("--packages", type=str, help="Comma-separated package list to check (overrides built-in list)")
    args = p.parse_args()

    python_exe = sys.executable
    if args.packages:
        packages = [x.strip() for x in args.packages.split(",") if x.strip()]
    else:
        packages = DEFAULT_PACKAGES

    missing = []
    for pkg in packages:
        ok = pip_show(python_exe, pkg)
        if not ok:
            missing.append(pkg)

    if not missing:
        print("All requested packages are present.")
    else:
        print("Missing packages:", ", ".join(missing))
        if args.install:
            failed = []
            for pkg in missing:
                print(f"Installing {pkg}...")
                ok = pip_install(python_exe, pkg)
                if not ok:
                    failed.append(pkg)
            if failed:
                print("Failed to install:", ", ".join(failed))
                sys.exit(1)
            else:
                print("All missing packages installed successfully.")
        else:
            print("Run with --install to attempt automatic installation.")
            sys.exit(1)

    if args.build_optsne:
        ok = build_optsne(python_exe)
        if not ok:
            sys.exit(3)

    if args.precompile:
        ok = run_precompile(python_exe)
        if not ok:
            sys.exit(2)


if __name__ == "__main__":
    main()
