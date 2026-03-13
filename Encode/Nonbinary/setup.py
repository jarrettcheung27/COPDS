from pathlib import Path

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ROOT_DIR = Path(__file__).resolve().parent

sources = [
    "fftqspa_bindings.cpp",
    "bch_codec.cpp",
    "BCJRQSPA.cpp",
    "QaryLDPC.cpp",
    "Mapper.cpp",
    "Random.cpp",
    "Interleaver.cpp",
    "FiniteField2.cpp",
    "util.cpp",
    "Qary_Gauss_E.cpp",
    "lsfr.cpp",
]

ext_modules = [
    Pybind11Extension(
        "fftqspa",
        sources=sources,
        include_dirs=[".", str(ROOT_DIR)],
        cxx_std=11,
        extra_compile_args=["/EHsc", "/D_CRT_SECURE_NO_WARNINGS", "/DNOMINMAX"],
    )
]

setup(
    name="fftqspa",
    version="0.1.0",
    description="FFTQSPA BiBo encoder/decoder bindings",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
