from Bio import SeqIO
from Bio.Seq import Seq
import string
from os import path
import pandas as pd

def load_seq(seq_spec):
    def isnucalphabet(seq_string):
        seq_string = seq_string.upper()
        exclude_set = list(string.ascii_uppercase) 
        exclude_set.remove('A')
        exclude_set.remove('C')
        exclude_set.remove('G')
        exclude_set.remove('T')
        exclude_set.remove('U')
        return not 1 in [c in seq_string for c in exclude_set]
    #check whether sec_spec could be the actual sequence
    #1. is it a BioPython sequence object?
    if type(seq_spec) == Seq:
        return str(seq_spec)
    elif type(seq_spec) == str:
        #2. is it a DNA sequence string?
        if len(seq_spec) > 60 and isnucalphabet(seq_spec):
            return seq_spec
        elif path.exists(seq_spec):
            #try parsing from one of the file types that SeqIO can handle
            filetypes = ['gb','fasta']
            for filetype in filetypes:
                parsed_records = list(SeqIO.parse(seq_spec, filetype))
                if len(parsed_records) > 0:
                      return str(parsed_records[0].seq)
            #try loading as simple .txt file
            with open(test, 'r') as reader:
                text_str = reader.read()
            if isnucalphabet(text_str):
                return text_str
    return

def extend_oligo(seq,start_idx,target_Tm=74):
        oli_seq = ''
        Tm = 0
        idx = start_idx
        while Tm < target_Tm:
            this_nt = seq[idx]
            if this_nt in ['A','T']:
                Tm += 2
            else:
                Tm += 4
            oli_seq += this_nt
            idx += 1
        return oli_seq
    
def rev_comp(seq,alphabet ='DNA'):
    comp_dict = {'A':'T','C':'G','G':'C','T':'A','U':'A'}
    if alphabet == 'RNA': 
        comp_dict['A'] = 'U'
    return ''.join([comp_dict[seq[n]] for n in range(len(seq)-1,-1,-1)])
    
def gibson_make_primers(vector_spec,target_spec, target_start_offset = 0,target_end_offset = 0,REs=[], homology_length=25, output='list'):
        
    """
    Generates a pair of sequences for primers that can be used to
    clone the target_file sequence into a vector_file sequence that has
    been linearised with one or two restriction enzymes using Gibson
    assembly.
    
    Parameters
    ----------
    vector_file: a genbank file containing the vector sequence
    target_file: a genbank file containing the sequence to be cloned into
                the vector sequence
    REs: list of valid names for one or two restriction enzymes used to
                linearize the vector sequence
    target_start_offset: the numbers of nucleotides to be excluded from
                the beginning of the target file (default = 0)
    target_end_offset: the numbers of nucleotides to be excluded from
                the end of the target file (default = 0)
    homology_length: the number of nucleotides to be added to the oligos
                that form the homology region for the Gibson assembly
                (default = 25).
    output: determines how the output is formatted. Allowed values 
                are "list" (default) and "fasta".
    
    Returns
    -------
    Oligonucleotide sequences in the specified format.
    """
    
    #check that output has an allowed value - set to default if not allowed
    if output not in ['fasta','list']:
        output = 'list'
    
    #read in vector sequence
    vector_seq = load_seq(vector_spec)
    
    #read in target sequence
    target_seq = load_seq(target_spec)
    
    #prepare package data for use
    try:
        import importlib.resources as pkg_resources
    except ImportError:
        # Try backported to PY<37 `importlib_resources`.
        import importlib_resources as pkg_resources
    from . import Data  # relative-import the *package* containing the data
    #import the stored data for S cerevisiae
    with open(Data.__path__[0] + '/RE_List.csv') as read_file:
        RE_ref = pd.read_csv(read_file)
    RE_ref.Name = RE_ref.Name.str.upper()
        
    #test whether correct number of restriction enzyme sites has been specified
    if type(REs) == list:
        REs = [RE.upper() for RE in REs]
        if len(REs) not in [1,2]:
            print('Only one or two restriction enzyme sites allowed.')
            return
        if len(REs) == 1:
            REs = REs + REs
        for RE in REs:
            if RE not in RE_ref.Name.values:
                print(RE + ' is an unknown restriction enzyme.')
                return
    else:
        if REs.upper() not in RE_ref.Name.values:
            print('Restriction enzyme name not recognised.')
            return
        REs = [REs.upper(),REs.upper()]
    
    #Import RE recognition sequences
    RE_sites = []
    for RE in REs:
        RE_sites.append(RE_ref.loc[RE_ref['Name'] == RE]['Motif'].values[0])
    RE_locations = []
    #if RE sequences exist exactly once in vector, record their index
    for RE_site in RE_sites: 
        if RE_site not in  vector_seq:
            print('Restriction sites for ' + RE + ' not found in vector')
            return
        if vector_seq.count(RE_site) > 1:
            print('Restriction site for ' + RE + ' found more than once in vector')
            return           
        RE_locations.append(vector_seq.index(RE_site))
    #assemble the forward oligo sequence
    homol_start = RE_locations[0]-homology_length
    homol_end = RE_locations[0]
    if homol_start < 0 :
        homology1 = vector_seq[homol_start:] + vector_seq[:homol_end]
    else:
        homology1 = vector_seq[homol_start:homol_end]
    target_anneal1 = extend_oligo(target_seq,target_start_offset)
    oligo1 = homology1 + RE_sites[0] + target_anneal1    
    #reverse complement the vector in order to process the reverse oligo correctly
    vector_seq_rc  = str(Seq(vector_seq).reverse_complement())
    target_seq_rc = str(Seq(target_seq).reverse_complement())
    #assemble the reverse oligo sequence
    RE_location_rc = vector_seq_rc.index(RE_sites[1])
    rc_homol_start = RE_location_rc-homology_length
    rc_homol_end = RE_location_rc
    if rc_homol_start < 0 :
        homology2 = vector_seq_rc[rc_homol_start:] + vector_seq_rc[:rc_homol_end]
    else:
        homology2 = vector_seq[rc_homol_start:rc_homol_end]
    target_anneal2 = extend_oligo(target_seq_rc,target_end_offset)
    oligo2 = homology2 + RE_sites[1] + target_anneal2
    
    #assemble the output
    if output=='fasta':
        return '>forward\n' + oligo1 + '\n>reverse\n' + oligo2
    elif output=='list':
        return[oligo1,oligo2]