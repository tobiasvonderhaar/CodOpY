from CodOpY.analyse import (plot_opt,translate,time_seq,report_repeats)
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
    
def test_translate():
    #define a DNA sequence containng all 64 codons
    seq = 'AAAAACAATAAGACAACCACTACGATAATCATTATGAGAAGCAGTAGGCAACACCATCAGCCACCCCCTCCGCTACTCCTTCTGCGACGCCGTCGGTAATACTATTAGTCATCCTCTTCGTTATTCTTTTTGTGATGCTGTTGGGAAGACGATGAGGCAGCCGCTGCGGTAGTCGTTGTGGGAGGCGGTGGG'
    #assert that thie translates into the correct amino acid sequence
    assert translate(seq) == 'KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG'
    
def test_time_seq():
    #define test sequence
    seq = 'GCGTGCGACGAATTTGGTCATATTAAATTGATGAACCCGCAAAGATCCACTGTCTGGTAC'
    #assert that time_seq returns the correct dictionary for this sequence for yeast
    assert time_seq(seq) == {'Decoding time': 2.2536, 'Average decoding time per codon': 0.11268, 'CV': 46.068781688038705}
    
def test_report_repeats():
    #define test sequence
    seq = 'MKRSAAAAATGCELLLLLLLLLLRRSTNNNQQQQ'
    #assert that the number of repeats reported for this sequence is correct
    assert report_repeats(seq) == {'A': (5, '@5'), 'L': (10, '@14'), 'Q': (4, '@31')}