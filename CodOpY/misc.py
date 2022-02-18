import pandas as pd

def initiate_code_table():
    '''
    This function returns a dataframe describing the standard genetic
    code.
    
    Parameters
    ==========
    The function does not accept any parameters.
    
    Returns
    =======
    pandas.core.frame.DataFrame
        A dataframe specifying the amino acids in one letter and three
        letter ode and their corresponding codons.
    
    '''

    codons = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC',
              'AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT',
              'CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC',
              'CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT',
              'GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC',
              'TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT', 
              'TTA','TTC','TTG','TTT']
    three_letter = ['Lys','Asn','Lys','Asn','Thr','Thr','Thr','Thr','Arg',
                    'Ser','Arg','Ser','Ile','Ile','Met','Ile','Gln','His',
                    'Gln','His','Pro','Pro','Pro','Pro','Arg','Arg','Arg',
                    'Arg','Leu','Leu','Leu','Leu','Glu','Asp','Glu','Asp',
                    'Ala','Ala','Ala','Ala','Gly','Gly','Gly','Gly','Val',
                    'Val','Val','Val','Stop','Tyr','Stop','Tyr','Ser','Ser',
                    'Ser','Ser','Stop','Cys','Trp','Cys','Leu','Phe','Leu',
                    'Phe']
    one_letter = ['K','N','K','N','T','T','T','T','R','S','R','S','I','I',
                  'M','I','Q','H','Q','H','P','P','P','P','R','R','R','R',
                  'L','L','L','L','E','D','E','D','A','A','A','A','G','G',
                  'G','G','V','V','V','V','*','Y','*','Y','S','S','S','S',
                  '*','C','W','C','L','F','L','F']
    return pd.DataFrame({'codon':codons,'three.letter':three_letter,
                         'one.letter':one_letter})

#=================================================================

def save_fasta(seqs, ids=[], filename='save_fasta.fa'):
    '''
    Saves an input sequence or set of input sequences in a specified 
    file in fasta format.
    
    Parameters
    ==========
    
    seqs: str or list of str
        the sequence or sequences to be saved.
    ids: str, int or float or list of these types
        the id(s) for the sequence(s) in seqs. If the length of ids 
        does not equal the length of seqs ids are ignored.
    filname: str
        the filename of the output file.
        
    Returns
    =======
    
    The function saves a file containing the sequences in seqs in 
    fasta format.
    
    '''
    if type(ids) != list:
        ids =[str(ids)]
    if type(seqs) == str:
            seqs = [seqs]
    if (ids == []) or (len(seqs) != len(ids)):
        ids = ['>'] * len(seqs)
    add_str = ''
    for idx,seq in enumerate(seqs):
        add_str += '>' + str(ids[idx]).replace('>','') + '\n'
        try:
            add_str += seq + '\n'
        except:
            raise ValueError('Invalid data type for sequence')
    with open(filename, 'w') as f:
        f.write(add_str)
        
#=================================================================

def load_from_Data(ref_table):
    
    '''
    This function loads auxiliary datasets provided with CodOpY in
    its 'Data' module.
    
    Parameters
    ==========
    
    ref_table : str
        the name of the filename to be loaded without the '.csv' 
        suffix.
        
    Returns
    =======
    
    pandas.core.frame.DataFrame
        the .csv data loaded into a pandas Dataframe.
    '''
    #prepare package data for use
    try:
        import importlib.resources as pkg_resources
    except ImportError:
        # Try backported to PY<37 `importlib_resources`.
        import importlib_resources as pkg_resources
    from CodOpY import Data  # import the *package* containing the data

    #import the stored data for the dataset in question
    try:
        with open(list(Data.__path__)[0] + '/' + ref_table + '.csv') as \
        read_file:
            frame = pd.read_csv(read_file)
    except:
        raise ValueError('ref_table does not refer to a valid dataset.')
    return frame

#=================================================================

def codon_choices(aa,ref_table='Scer',parameter='decoding.time'):
    '''
    Reports the codon choices for a specified amino acid.

    Parameters:
    ===========
    aa : str
        An amino acid in one or three letter code.
    ref_table : str
        A valid name for a parameterset such as 'Scer'.
    parameter : str
        The name of the parameter to be displayed alongside the codons, which must
        correspond to one of the columns of the specified parset.

    Returns:
    ========
        pandas.core.frame.DataFrame
    '''
    ref_table = load_from_Data(ref_table)
    if len(aa) == 1:
        aa = aa.upper()
        if aa not in ref_table['one.letter'].values:
            raise ValueError('Invalid amino acid choice.')
        select_table = ref_table.loc[ref_table['one.letter'] == aa]
    elif len(aa) == 3:
        aa = aa[0].upper() + aa[1:].lower()
        if aa not in ref_table['three.letter'].values:
            raise ValueError('Invalid amino acid choice.')
        select_table = ref_table.loc[ref_table['three.letter'] == aa]
    else:
        raise ValueError('Invalid amino acid choice.')
    return select_table[['codon','three.letter','one.letter',parameter]].sort_values(by=parameter).reset_index(drop=True)