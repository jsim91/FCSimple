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
  python fcs_install_python_deps.py [--install] [--precompile] [--packages pkg1,pkg2]

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
    "openTSNE"
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

    if args.precompile:
        ok = run_precompile(python_exe)
        if not ok:
            sys.exit(2)


if __name__ == "__main__":
    main()
