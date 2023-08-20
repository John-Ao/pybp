from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.1.0"

setup(
    name="pybp",
    version=__version__,
    author="aowenxuan",
    author_email="aowenxuan@tsingroc.com",
    url="https://git.tsingroc.com/or/pybp",
    description="2D Bin-Packing",
    long_description="",
    ext_modules=[
        Pybind11Extension(
            "pybp", ["main.cc"], define_macros=[("VERSION_INFO", __version__)]
        )
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
    include_package_data=True,
)
