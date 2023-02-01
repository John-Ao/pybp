"""
    2D Bin Packing
    -----------------------

    .. currentmodule:: pybp

    .. autosummary::
        :toctree: _generate

        get_bin_out_fast_1
        get_bin_out_fast_2
    """
from __future__ import annotations
import pybp
import typing
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "bp_max_rect",
    "get_bin_out_fast_1",
    "get_bin_out_fast_2",
    "get_bin_out_fast_3",
    "get_bin_out_fast_4"
]


def bp_max_rect(items: numpy.ndarray[numpy.float64], rotatable: bool, bin_width: float = 2440, bin_height: float = 1220) -> typing.Tuple[typing.List[int], typing.List[float]]:
    """
    Do MaxRect bin-packing algorithm on the given input

    items: float np.array, Nx2

    rotatable: bool, whether the items are rotatable (IGNORED)

    bin_width: float, width of the bins

    bin_height: float, height of the bins

    return: (bin_id of each item, usage of each bin)
    """
def get_bin_out_fast_1(items: numpy.ndarray[numpy.float64], rotatable: bool, bin_width: float = 2440, bin_height: float = 1220) -> int:
    """
    Do bin-packing algorithm on the given input

    items: float np.array, Nx2

    rotatable: bool, whether the items are rotatable

    bin_width: float, width of the bins

    bin_height: float, height of the bins

    return: # of used bins
    """
def get_bin_out_fast_2(items: numpy.ndarray[numpy.float64], rotatable: bool, bin_width: float = 2440, bin_height: float = 1220) -> typing.Tuple[int, float]:
    """
    Do bin-packing algorithm on the given input

    items: float np.array, Nx2

    rotatable: bool, whether the items are rotatable

    bin_width: float, width of the bins

    bin_height: float, height of the bins

    return: (# of used bins, % of usage of the last bin)
    """
def get_bin_out_fast_3(items: numpy.ndarray[numpy.float64], rotatable: bool, bin_width: float = 2440, bin_height: float = 1220) -> typing.Tuple[typing.List[int], typing.List[float]]:
    """
    Do bin-packing algorithm on the given input

    items: float np.array, Nx2

    rotatable: bool, whether the items are rotatable

    bin_width: float, width of the bins

    bin_height: float, height of the bins

    return: (bin_id of each item, usage of each bin)
    """
def get_bin_out_fast_4(items: numpy.ndarray[numpy.float64], rotatable: bool, bin_width: float = 2440, bin_height: float = 1220) -> typing.Tuple[typing.List[int], typing.List[float]]:
    """
    Do bin-packing algorithm on the given input

    items: float np.array, Nx2

    rotatable: bool, whether the items are rotatable

    bin_width: float, width of the bins

    bin_height: float, height of the bins

    return: (bin_id of each item, usage of each bin)
    """
__version__ = '0.0.8'
