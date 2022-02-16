from CodOpY.plot import (plot_opt)
import matplotlib
import pytest

def test_plot_opt():
    #prepare a test figure
    seq = "ATGGATCGATCGATGCAGTCGATCGTAGCTGACTCATGCTACATGGATCGATCGATGCAGTCGATCGTAGCTGACTTATGCTAC"
    test = plot_opt(seq)
    #assert that the function returns a tuple of maplotlib figure and axes
    assert isinstance(test[0],matplotlib.figure.Figure)
    #assert that labels and axis limits are as predicted
    assert str(test[0].get_axes()) == "[<AxesSubplot:xlabel='Codon No', ylabel='decoding.time'>]"
    assert test[1].get_xlim() == pytest.approx((-1.35, 28.35))