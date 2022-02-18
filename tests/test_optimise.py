from CodOpY.optimise import (opt_seq, remove_RE,Parset)



def test_opt_seq():
    #define a string containing all 20 single letter amino acid codes
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    #assert that this produces the correct optimised sequence for yeast
    assert opt_seq(AA,diversify=[]) == 'GCGTGCGACGAATTTGGTCATATTAAATTGATGAACCCGCAAAGATCCACTGTCTGGTAC'
    #assert that when used for yeast with the default diversify parameters this never includes slow codons
    for rep in range(25): 
        assert 'GGG' not in opt_seq(AA)
        assert 'GGA' not in opt_seq(AA)
    #test the enforce option
    assert opt_seq('MKKKGTACK',enforce={'K':'AAC'}) == 'ATGAACAACAACGGTACTGCGTGCAAC'

def test_remove_RE():
   #define test sequence
    mcs ='AATTCCATATGTTAATTAAGGCGCGCCCAATTGGATCCA'
   #assert that this returns the correct sequence where RE has been removed
    assert remove_RE('NdeI',mcs) == 'AATTCCATTTGTTAATTAAGGCGCGCCCAATTGGATCCA'

def test_Parset():
    parset = Parset()
    assert parset.__repr__() == 'This is a Parset instance containing a CodOpY parameterset.'
	