import os
import subprocess
import sys
from pathlib import Path

from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext


def read_requirements(filename):
    with open(filename) as req_file:
        return req_file.read().splitlines()


class CMakeBuildExt(build_ext):
    """Custom build_ext that compiles:
    1. The janas_core Python C extension (via setuptools)
    2. The standalone C++ apps (janas_app_meanMinMax, janas_app_starProcess) via CMake
    """

    def run(self):
        target_dir = os.path.join(self.build_lib, "janas")
        os.makedirs(target_dir, exist_ok=True)
        super().run()
        self._build_cmake_apps()
        # For editable/inplace installs, also copy to source tree
        if self.inplace:
            src_bin = Path("src") / "janas" / "bin"
            lib_bin = Path(self.build_lib) / "janas" / "bin"
            if lib_bin.exists():
                src_bin.mkdir(parents=True, exist_ok=True)
                for binary in lib_bin.iterdir():
                    if binary.is_file():
                        dest = src_bin / binary.name
                        dest.write_bytes(binary.read_bytes())
                        dest.chmod(0o755)

    def copy_extensions_to_source(self):
        target_dir = os.path.join("src", "src", "janas")
        os.makedirs(target_dir, exist_ok=True)
        super().copy_extensions_to_source()
        # Also copy compiled C++ apps to source tree for editable installs
        src_bin = Path("src") / "janas" / "bin"
        lib_bin = Path(self.build_lib) / "janas" / "bin"
        if lib_bin.exists():
            src_bin.mkdir(parents=True, exist_ok=True)
            for binary in lib_bin.iterdir():
                if binary.is_file():
                    dest = src_bin / binary.name
                    dest.write_bytes(binary.read_bytes())
                    dest.chmod(0o755)

    def _build_cmake_apps(self):
        source_dir = Path(__file__).resolve().parent
        build_dir = Path(self.build_temp) / "cmake_build"
        build_dir.mkdir(parents=True, exist_ok=True)

        bin_dir = Path(self.build_lib) / "janas" / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_INSTALL_PREFIX={bin_dir}",
            f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={bin_dir}",
            "-DCMAKE_BUILD_TYPE=Release",
        ]

        subprocess.check_call(
            ["cmake", str(source_dir)] + cmake_args,
            cwd=str(build_dir),
        )
        subprocess.check_call(
            ["cmake", "--build", ".", "--parallel"],
            cwd=str(build_dir),
        )

        # Copy binaries to the package bin/ directory
        app_bin = build_dir / "app_bin"
        if app_bin.exists():
            for binary in app_bin.iterdir():
                if binary.is_file() and os.access(str(binary), os.X_OK):
                    dest = bin_dir / binary.name
                    dest.write_bytes(binary.read_bytes())
                    dest.chmod(0o755)


# Python C extension
janas_core = Extension(
    'janas.janas_core',
    sources=['src/src_cpp/janas_core.cpp'],
    language='c++',
)

setup(
    name='janas',
    version='2.1.5',
    description='Pipeline for per-particle selection and 3D class reassignment in single-particle cryo-EM.',
    author='Mauro Maiorca',
    author_email='mauro.maiorca@cssb-hamburg.de',
    url='https://gitlab.com/topf-lab/janas',
    python_requires='>=3.8',
    install_requires=read_requirements('requirements.txt'),
    extras_require={
        # GPU support for janas_reconstructor (used when --noExternalPrograms
        # --gpu N is set). Pulls torch from the default PyPI index. For a
        # specific CUDA build, install torch manually with --index-url
        # https://download.pytorch.org/whl/cuXYZ before/after this package.
        'gpu': ['torch>=2.0'],
    },
    packages=find_packages(where='src', include=['janas*']),
    package_dir={'': 'src'},
    ext_modules=[janas_core],
    cmdclass={'build_ext': CMakeBuildExt},
    include_package_data=True,
    package_data={'janas': ['config.json', 'bin/*', 'images/*.png']},
    entry_points={
        'console_scripts': [
            'janas = janas.janas_cmd_caller:main',
            'janas_utils = janas.janas_cmd_utils:main',
            'janas_optimizer = janas.janas_cmd_optimizer:main',
            'janas_session_manager = janas.janas_cmd_session_manager:main',
            'janas_reconstructor = janas.janas_reconstructor:main',
            'janas_app_meanMinMax = janas.bin_wrapper:meanMinMax_main',
            'janas_app_starProcess = janas.bin_wrapper:starProcess_main',
        ],
    },
    license='MIT',
    zip_safe=False,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
