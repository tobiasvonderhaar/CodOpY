import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from CodOpY.optimise import opt_seq
from CodOpY.misc import load_from_Data



def translate(seq):
    '''
    Returns an amino acid sequence translated for an input DNA 
    sequence.

    Parameters
    ==========
    seq : str
        The DNA sequence to be translated.

    Returns
    str
        The amino acid sequence translated from seq.
    '''

    codedict = {'AAA':'K','AAC':'N','AAG':'K','AAT':'N','ACA':'T','ACC':'T',
                'ACG':'T','ACT':'T','AGA':'R','AGC':'S','AGG':'R','AGT':'S',
                'ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H',
                'CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
                'CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L',
                'CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E','GAT':'D',
                'GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G',
                'GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V','GTT':'V',
                'TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S',
                'TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W','TGT':'C',
                'TTA':'L','TTC':'F','TTG':'L','TTT':'F'}

    if type(seq) == str and all([n in ['A','C','T','G'] for n in seq.upper()]):
        seq = [seq[i:i+3] for i in range(0,len(seq),3) if len(seq[i:i+3]) == 3]
    elif type(seq) == list and all([len(c) == 3 for c in seq]):
        pass
    else:
        raise ValueError('This sequence is not recognised as a valid DNA or' +
                          ' codon sequence.')

    return ''.join([codedict[c] for c in seq])

#================================================================================

def time_seq(input_seq,ref_table='Scer'):
    '''time_seq calculates the time it takes to decode a DNA or RNA sequence.

    Parameters
    ----------
    input_seq : str
        The DNA or RNA sequence for which the tiem properties are beiing returned.
    ref_table : str
        The name of the reference table to be used, eg 'Scer' for S. cerevisiae

    Returns
    -------
    dict
        A dictionary containing the overall Decoding time in seconds,
        the Average decoding time per codon in seconds, and the CV
        (coefficient of variation of the decoding time per codon)
    '''

    #convert the input sequence to DNA if RNA
    input_seq = input_seq.replace('U','T')
    #make a look-up dictionary of the decoding times available for an amino acid
    #from the reference table
    ref = load_from_Data(ref_table)
    ref['codon'] = ref['codon'].str.replace('U','T')
    time_dict_by_codon = dict(zip(ref['codon'], ref['decoding.time']))
    #add pseudo data for stop codons
    time_dict_by_codon['TGA'] = 0
    time_dict_by_codon['TAA'] = 0
    time_dict_by_codon['TAG'] = 0
    #calculate the decoding time properties of the sequence
    codon_seq = [input_seq[n:n+3] for n in range(0,len(input_seq),3)]
    times_vec = [time_dict_by_codon[codon] for codon in codon_seq]
    results = {}
    results['Decoding time'] = np.sum(times_vec)
    results['Average decoding time per codon'] = results['Decoding time'] / len(codon_seq)
    results['CV'] = np.std(times_vec, ddof=1) / np.mean(times_vec) * 100
    return results

#================================================================================

def report_repeats(seq, threshold=4):
    '''
    Reports the repeat length for any single amino acid repeat longer than the threshold value.

    Parameters
    ==========
    seq : str
        The amino acid sequence to be analysed
    threshold : int
        The minimum number of consecutive amino acids reported as a repeat. Default = 4.

    Returns
    =======
    
    dict
        A dictionary of amino acids for which repeats were found and the longest repeat length.
        Returns an empty dictionary if no repeats were found.
    '''

    ret_dict = {}
    for aa in seq:
        l = threshold
        while aa * l in seq:
            ret_dict[aa] = (l,'@')
            l += 1
    for key,value in ret_dict.items():
        ret_dict[key] = (value[0],value[1] + str(seq.index(key*value[0]) + 1))
    return ret_dict

#================================================================================

def test_RE(RE, test_seq):
    '''
    Tests whether an RE site is present in a sequence.

    Parameters
    ==========
    REs : str
        The name of the enzyme, or the sequence, to be tested.

    test_seq : str
        The DNA sequence to be tested.

    Returns
    =======
    list
        A list of names of those enzymes for which sites were found in
        test_seq, or an empty list of none of the enzyme sites was found.
    '''

    RE_ref = load_from_Data('RE_List')
    
    #convert RE names to upper case for comparison to 'RE' variable
    RE_ref.Name = RE_ref.Name.str.upper()

    #is site a restriction enzyme name?
    RE = RE.upper()
    if RE in RE_ref.Name.values:
        RE_seq = RE_ref.loc[RE_ref['Name'] == RE]['Motif'].values[0]
    #else is site a valid DNA sequence?
    elif not 0 in [c in ['A','C','T','G','W','S','M','K','R','Y','N'] for c in RE]:
        RE_seq = RE
    else:
        print('No known Restriction enzyme site or valid DNA sequence specified.')
        return
    #remove leading or trailing Ns from RE site
    while RE_seq[0] == 'N':
        RE_seq = RE_site[1:]
    while RE_seq[-1] == 'N':
        RE_seq = RE_site[:-1]
    #does RE_site now have more than 5 'N' (this becomes very inefficient)
    if RE_seq.count('N') > 5:
        print('Too many N - the maximum number of N allowed in the RE sequence is 5')
        return

    #does the RE site contain ambiguous nucleotide symbols?
    #If so list all corresponding unambiguous sequences
    ambiguous_bases = {'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],'R':['A','G'],'Y':['C','T'],
                       'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}
    unamb_RE_seqs = [RE_seq]
    while 0 in [c in ['A','C','T','G'] for c in unamb_RE_seqs[0]]:
        new_options = []
        for seq in unamb_RE_seqs:
            for idx, nt in enumerate(seq):
                if nt in ambiguous_bases.keys():
                    for nt in ambiguous_bases[nt]:
                        new_options.append(seq[:idx] + nt + seq[idx+1:])
                    break
        unamb_RE_seqs = new_options


    #convert sequence to valid uppercase DNA sequence
    test_seq = test_seq.upper().replace('U','T')

    #is site in seq?
    if not 1 in [R in test_seq for R in unamb_RE_seqs]:
        print(RE + ' not found in this sequence. \n')
        return
    #if yes go through the sequence and record all sites of occurrence
    else:
        found_sites = []
        sub_seq = test_seq
        while 1 in [R in sub_seq for R in unamb_RE_seqs]:
            for R in unamb_RE_seqs:
                if R in sub_seq:
                    if len(found_sites) >= 1:
                        index_shift = found_sites[-1]
                    else:
                        index_shift = 0
                    found_sites.append(sub_seq.find(R)+ index_shift)
                    sub_seq = sub_seq[found_sites[-1]+1:]
        print(RE + ' found at the following site(s): ' + str(found_sites).replace('[','').replace(']','') + '\n')
        return
#================================================================================

def count_aas(seq):
    '''Produces a bar graph of amino acid counts in a protein sequence'''

    aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    counts = []
    for aa in aas:
        counts.append(seq.count(aa))
    fig,ax = plt.subplots()
    ax.bar(aas,counts)
    ax.set_ylabel('counts')
    ax.set_xlabel('amino acid')
    return fig,ax

#==================================================================================================

def plot_opt(seq, ref_table = 'Scer', plot_par = 'decoding.time',window=25,plot_colors = ['gold','cornflowerblue','mediumpurple']):

    #convert sequence to DNA
    seq = seq.upper()
    seq = seq.replace('U','T')

    parameterset = load_from_Data(ref_table)

    #generate a lookup dictionary for the plotted parameter
    pardict = {}
    for c in parameterset['codon']:
        this_par = parameterset.loc[parameterset['codon'] == c][plot_par].values[0]
        this_c= c.replace('U','T')
        pardict[this_c] = this_par
    #make the extreme arameter sequences
    aaseq = translate(seq)
    max_seq = opt_seq(aaseq,ref_table = ref_table,optimise_by = [plot_par,max],diversify = [])
    min_seq = opt_seq(aaseq,ref_table = ref_table,optimise_by = [plot_par,min], diversify = [])

    codon_seq = [seq[i:i+3] for i in range(0,len(seq),3) if len(seq[i:i+3]) == 3]
    max_codon_seq = [max_seq[i:i+3] for i in range(0,len(max_seq),3) if len(max_seq[i:i+3]) == 3]
    min_codon_seq = [min_seq[i:i+3] for i in range(0,len(min_seq),3) if len(min_seq[i:i+3]) == 3]


    #prepare vectors of the codon parameter for each sequence:
    seq_pars = [pardict[c] for c in codon_seq]
    max_pars = [pardict[c] for c in max_codon_seq]
    min_pars = [pardict[c] for c in min_codon_seq]
    x = list(range(len(seq_pars)))

    smooth_seq_pars = [sum(seq_pars[i:i+window])/window for i in range(len(seq_pars) - window)]
    smooth_max_pars = [sum(max_pars[i:i+window])/window for i in range(len(max_pars) - window)]
    smooth_min_pars = [sum(min_pars[i:i+window])/window for i in range(len(min_pars) - window)]

    smooth_x_offset = int(window/2)
    smooth_x = list(range(smooth_x_offset,smooth_x_offset + len(smooth_min_pars)))
    #prepare the plot
    fig,ax = plt.subplots()
    ax.scatter(x, seq_pars,c=plot_colors[0],s=5,alpha=0.5)
    ax.plot(smooth_x,smooth_seq_pars,c=plot_colors[0],label='actual')
    ax.plot(smooth_x,smooth_max_pars,c=plot_colors[1],label='max')
    ax.plot(smooth_x,smooth_min_pars,c=plot_colors[2],label='min')

    ax.set_xlabel('Codon No')
    ax.set_ylabel(plot_par)

    plt.legend(loc='upper right')

    return fig,ax
