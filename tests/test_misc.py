from CodOpY.misc import (initiate_code_table, save_fasta)

import pandas
import os

def test_initiate_code_table():
    v = initiate_code_table()
    assert list(v.columns) == ['codon', 'three.letter', 'one.letter']
    assert v.loc[v['codon'] == 'ATG']['one.letter'].values[0] == 'M'
    assert v.loc[v['codon'] == 'TGG']['one.letter'].values[0] == 'W'
    assert v.loc[v['codon'] == 'CAG']['one.letter'].values[0] == 'Q'
    assert v.loc[v['codon'] == 'GGG']['one.letter'].values[0] == 'G'
    
def test_save_fasta():
    filename = 'tst.txt'
    save_fasta(['ATGATGAT','TAGCCGCGATC'],ids=[1,'>test'], filename=filename)
    with open(filename, 'r') as f:
        re_read = f.read()
    os.remove(filename)
    assert re_read == '>1\nATGATGAT\n>test\nTAGCCGCGATC\n'
    
    save_fasta('ATGATGAT',ids='test', filename=filename)
    with open(filename, 'r') as f:
        re_read = f.read()
    os.remove(filename)
    assert re_read == '>test\nATGATGAT\n'