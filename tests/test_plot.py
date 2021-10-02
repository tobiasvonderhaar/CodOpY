from CodOpY.plot import (plot_opt)

def test_plot_opt():
    #prepare a test figure
    seq = "ATGGATCGATCGATGCAGTCGATCGTAGCTGACTCATGCTACATGGATCGATCGATGCAGTCGATCGTAGCTGACTTATGCTAC"
    test = plot_opt(seq)
    #assert that the function returns maplotlib figure and axes
    assert str(type(test[0])) == "<class 'matplotlib.figure.Figure'>"
    assert str(type(test[1])) == "<class 'matplotlib.axes._subplots.AxesSubplot'>"
    #assert that labels and axis limits are as predicted
    assert str(test[0].get_axes()) == "[<AxesSubplot:xlabel='Codon No', ylabel='decoding.time'>]"
    assert test[1].get_xlim() == (-1.35, 28.35)