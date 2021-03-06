from CodOpY.optimise import (translate, opt_seq, remove_RE,time_seq,retrieve_kazusa,report_repeats)

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
    #assert that this produces the correct optiise sequence for HEK293 cells
    assert opt_seq(AA,diversify=[],ref_table='HEK') == 'GCTTGCGACGAATTCGGACACATTAAACTTATGAACCCTCAAAGATCTACTGTTTGGTAC'
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
    assert time_seq(seq,ref_table='HEK') == {'Decoding time': 23.310745777, 'Average decoding time per codon': 1.16553728885, 'CV': 33.04953199594406}

def test_retrieve_kazusa():
    #assert that the first value returned for the Kazua data for yeast is correct
    assert retrieve_kazusa(4932)['usage.frequency'][0] == 41.9

def test_report_repeats():
    #define test sequence
    seq = 'MKRSAAAAATGCELLLLLLLLLLRRSTNNNQQQQ'
    #assert that the number of repeats reported for this sequence is correct
    assert report_repeats(seq) == {'A': (5, '@5'), 'L': (10, '@14'), 'Q': (4, '@31')}
