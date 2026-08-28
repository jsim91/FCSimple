import os
import shutil
import subprocess
import sysconfig
from os import path

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

PACKAGE_NAME = "MulticoreTSNE"
VERSION = "0.1"


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        # Prefer an explicitly provided cmake on PATH or in the Python Scripts
        # directory (e.g. installed via `pip install cmake`).
        cmake = shutil.which("cmake")
        if cmake is None:
            scripts = sysconfig.get_path("scripts")
            cand = path.join(scripts, "cmake.exe") if os.name == "nt" else path.join(scripts, "cmake")
            if path.exists(cand):
                cmake = cand
        if cmake is None:
            raise RuntimeError("cmake not found. Install it, e.g. `pip install cmake`.")

        for ext in self.extensions:
            self._build_one(ext, cmake)

    def _build_one(self, ext, cmake):
        source_dir = ext.sourcedir
        build_temp = self.build_temp
        ext_path = self.get_ext_fullpath(ext.name)
        ext_dir = path.abspath(path.dirname(ext_path))

        shutil.rmtree(build_temp, ignore_errors=True)
        os.makedirs(build_temp, exist_ok=True)

        build_type = "Debug" if self.debug else "Release"
        # Set both generic and configuration-specific output dirs so the
        # compiled extension lands next to the Python package regardless of
        # generator/configuration.
        lib_dir = "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(ext_dir)

        configure = [
            cmake,
            "-S", source_dir,
            "-B", build_temp,
            "-DCMAKE_BUILD_TYPE={}".format(build_type),
            lib_dir,
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG={}".format(ext_dir),
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE={}".format(ext_dir),
        ]

        generator = os.environ.get("CMAKE_GENERATOR")
        if generator:
            configure += ["-G", generator]

        compiler = os.environ.get("CMAKE_CXX_COMPILER")
        if compiler:
            configure += ["-DCMAKE_CXX_COMPILER={}".format(compiler)]

        print("Running CMake configure:", " ".join(configure))
        subprocess.check_call(configure, cwd=build_temp)

        build_cmd = [cmake, "--build", build_temp]
        # MSVC requires the configuration be specified at build time; harmless
        # for other generators.
        if os.name == "nt":
            build_cmd += ["--config", build_type]
        print("Running CMake build:", " ".join(build_cmd))
        subprocess.check_call(build_cmd, cwd=build_temp)


if __name__ == "__main__":
    setup(
        name=PACKAGE_NAME,
        version=VERSION,
        description="Multicore version of t-SNE algorithm (opt-SNE fork).",
        author="Dmitry Ulyanov (based on L. Van der Maaten's code)",
        author_email="dmitry.ulyanov.msu@gmail.com",
        url="https://github.com/omiq-ai/Multicore-opt-SNE",
        install_requires=[
            "numpy",
            "cffi",
        ],
        packages=find_packages(),
        include_package_data=True,
        ext_modules=[CMakeExtension("MulticoreTSNE.MulticoreTSNE", sourcedir="multicore_tsne")],
        cmdclass={"build_ext": CMakeBuild},
    )