from CodOpY.optimise import (translate, opt_seq, remove_RE,time_seq,retrieve_kazusa,report_repeats,codon_choices)

def test_translate():
    #define a DNA sequence containng all 64 codons
    seq = 'AAAAACAATAAGACAACCACTACGATAATCATTATGAGAAGCAGTAGGCAACACCATCAGCCACCCCCTCCGCTACTCCTTCTGCGACGCCGTCGGTAATACTATTAGTCATCCTCTTCGTTATTCTTTTTGTGATGCTGTTGGGAAGACGATGAGGCAGCCGCTGCGGTAGTCGTTGTGGGAGGCGGTGGG'
    #assert that thie translates into the correct amino acid sequence
    assert translate(seq) == 'KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG'

def test_opt_seq():
    #define a string containing all 20 single letter amino acid codes
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    #assert that this produces the correct optimised sequence for yeast
    assert opt_seq(AA,diversify=[]) == 'GCGTGCGACGAATTTGGTCATATTAAATTGATGAACCCGCAAAGATCCACTGTCTGGTAC'
    #assert the when used for yeast with the default diversify parameters this never includes slow codons
    for rep in range(25):
        assert 'GGG' not in opt_seq(AA)
        assert 'GGA' not in opt_seq(AA)

def test_remove_RE():
   #define test sequence
    mcs ='AATTCCATATGTTAATTAAGGCGCGCCCAATTGGATCCA'
   #assert that this returns the correct sequence where RE has been removed
    assert remove_RE('NdeI',mcs) == 'AATTCCATTTGTTAATTAAGGCGCGCCCAATTGGATCCA'

def test_time_seq():
    #define test sequence
    seq = 'GCGTGCGACGAATTTGGTCATATTAAATTGATGAACCCGCAAAGATCCACTGTCTGGTAC'
    #assert that time_seq returns the correct dictionary for this sequence for yeast
    assert time_seq(seq) == {'Decoding time': 2.2536, 'Average decoding time per codon': 0.11268, 'CV': 46.068781688038705}

def test_retrieve_kazusa():
    #assert that the first value returned for the Kazusa data for yeast is correct
    assert retrieve_kazusa(4932)['usage.frequency'][0] == 41.9

def test_report_repeats():
    #define test sequence
    seq = 'MKRSAAAAATGCELLLLLLLLLLRRSTNNNQQQQ'
    #assert that the number of repeats reported for this sequence is correct
    assert report_repeats(seq) == {'A': (5, '@5'), 'L': (10, '@14'), 'Q': (4, '@31')}

def test_codon_choices():
    #tests are run with the default organism 'Scer'
    #assert that the function returns a pandas dataframe
    assert str(type(codon_choices("R"))) == "<class 'pandas.core.frame.DataFrame'>"
    #assert that the dataframe for test aino acids has the correct shape
    assert codon_choices("R").shape == (6,4)
    assert codon_choices("Met").shape == (1,4)
    codon_choices("Ile").shape == (3,4)
	