import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import requests
from CodOpY.misc import initiate_code_table, load_from_Data

class Parset():
    def __init__(self):
        self.ref_table = 'Scer'
        self.objective = 'decoding.time'
        self.objective_function = min
        self.diversify = []
        self.diversify_range = 2
        self.diversify_amino_acid_repeats = True
        self.enforce = {}
        self.exclude_REs = []
        self.reduce_speed_at = ()
        self.reduce_speed_by = 0
        
    def __repr__(self):
        return 'This is a Parset instance containing a CodOpY parameterset.'
    
    def __str__(self):
        printable_attributes = [item for item in dir(self) if '__' not in item]
        print('The parameters in this Parset instance are \n')
        for attr in printable_attributes:
            print(attr + ': ' + str(getattr(self,attr)))
        return ''

#==================================================================================================

def opt_seq(seq,diversify=[],
            diversify_range=0.2,ref_table = 'Scer',
            enforce={},
            optimise_by=['decoding.time',min]):

    '''
    Makes an optimised DNA sequence corresponding to an input amino
    acid sequence.

    Parameters
    ==========
    seq : str
        The amino aid sequence to be otimised

    diversify : list of str
        A list specifying individual amino acids for which codons
        should be diversified.

    diversify_range : float
        The proportion over which the optimisation parameters that is
        allowed to vary for diversified amno acids, relative to 1.

    ref_table : str
        the name of the data file containing the codon data.
        
    enforce : dict
        A dictionary for manually specifying codons for certain amino 
        acids (for these the optimisation parameters are overridden).
        Example: enforce = {'K':'AAG'} will always use AAG to encode 
        lysine (K).

    optimise_by : list of str and function
        The str part of optimise_by specifies which cloumn of the 
        ref_table shold be used for optimisation. Function can be min
        or max and specifies whether the highest or lowest value should
        be selected for optimisation.

    Returns
    =======
    str
        Returns a DNA sequence string.
    '''
    
    parameterset = load_from_Data(ref_table)

    #define the reverse translation dictionary: 
    #which codons should be considered for which amino acid?
    #if an amino acid is not in the diversify list, use the fastest codon
    #if an amino acid is in the diversify list, diversify codon choice by 
    #random draw from all codons for which the diversify_range threshold 
    #applies
    reverse_dict = {}
    for aa in parameterset['one.letter'].unique():
        codons_for_aa = parameterset.loc[parameterset['one.letter']==aa]
        if aa in diversify:
            if optimise_by[1] == min:
                diversify_threshold = min(codons_for_aa[optimise_by[0]]) * (1 + diversify_range)
                acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']<diversify_threshold]['codon'])
            elif optimise_by[1] == max:
                diversify_threshold = max(codons_for_aa[optimise_by[0]]) * (1 - diversify_range)
                acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']>diversify_threshold]['codon'])
            else:
                raise ValueError("Invalid amino acid specified in diversify")
        else:
            best_decoding_time_for_aa = optimise_by[1](codons_for_aa['decoding.time'])
            acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']==best_decoding_time_for_aa]['codon'])
        reverse_dict[aa] = acceptable_codons_for_aa
    #have codons been specified using the 'enforce' parameter?
    if len(enforce) > 0:
           for key, value in enforce.items():
                codon_set = list(parameterset.loc[parameterset['one.letter'] == key]['codon'].values)
                if value not in codon_set:
                    print('Non-standard genetic code enforced!')
                print('Codon for ' + key + ' enforced as ' + value)
                if type(value) != list:
                    value = [value]
                reverse_dict[key] = value
    codon_seq = []
    for seq_aa in seq:
        codon_seq.append(random.choice(reverse_dict[seq_aa]))
    return ''.join(codon_seq).replace('U','T')


#==================================================================================================

def remove_RE(site, test_seq, ref_table='Scer',optimise_by=['decoding.time',min],suppress_not_found=False):

    '''Removes restriction enzyme sites from DNA sequences without altering the encoded
    amino acid sequence and while maintaining codon optimisation as much as possible.

    Parameters
    ==========
    site : str
        the name of the restriction enzyme for which sites should be removed.

    test_seq : str
        the sequence from which sites are to be removed.

    ref_table : str
        the name of the reference table from which the optimisation information
        is being used.

    optimise_by : list of str and Function
        As for opt_seq, the name of the column of ref_table from which
        optimisation info is generated, and whether optimal is the minimum or
        maximum.

    Returns
    =======
    str
        A DNA sequence string.
    '''
    
    #load auxiliary package data
    RE_ref = load_from_Data('RE_List')
    codons = load_from_Data(ref_table)
    
    #convert RE names to upper case for comparison to 'RE' variable
    RE_ref.Name = RE_ref.Name.str.upper()

    #is site a restrictions enzyme name?
    site = site.upper()
    if site in RE_ref.Name.values:
        RE_seq = RE_ref.loc[RE_ref['Name'] == site]['Motif'].values[0]
    #is site a valid DNA sequence?
    elif not 0 in [c in ['A','C','T','G','W','S','M','K','R','Y','N'] for c in site]:
        RE_seq = site
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
        if not suppress_not_found:
            print(site + ' not found in this sequence')
        return
    new_seq = test_seq

    #prepare auxiliary data structures
    codon_pos = list(range(0,len(test_seq)-2,3))
    codons['codon'] = codons['codon'].str.replace('U','T')
    code_lookup = pd.Series(codons['one.letter'].values,index=codons.codon).to_dict()
    speed_lookup = pd.Series(codons[optimise_by[0]].values,index=codons.codon).to_dict()
    reverse_code_lookup = {}
    for aa in codons['one.letter'].unique():
        this_subset = codons.loc[codons['one.letter'] == aa]
        reverse_code_lookup[aa] = list(this_subset['codon'])

    #go through each RE site and replace one codon to remove it
    while 1 in [R in new_seq for R in unamb_RE_seqs]:
        #determine the starting nt of the first instance of the RE_site
        for this_RE_seq in unamb_RE_seqs:
            try:
                found_RE_index = new_seq.index(this_RE_seq)
                break
            except:
                pass
        #isolate the sequence of codons that contains the site
        subseq_start = max([x for x in codon_pos if found_RE_index >= x])
        if (found_RE_index + len(unamb_RE_seqs[0])) > len(new_seq)-2:
            subseq_stop = len(new_seq)
        else:
            subseq_stop = min([x for x in codon_pos if (found_RE_index + len(unamb_RE_seqs[0]) <= x)])
        subseq_contains_site= new_seq[subseq_start:subseq_stop]
        subseq_codons = [subseq_contains_site[n:n+3] for n in range(0,len(subseq_contains_site),3)]

        seq_vec,times_vec = [],[]
        for idx,sub_codon in enumerate(subseq_codons):
            this_aa = code_lookup[sub_codon]
            for alternative in reverse_code_lookup[this_aa]:
                if alternative != sub_codon:
                    new_subseq_codons = subseq_codons[:idx] + [alternative] + subseq_codons[idx+1:]
                    new_subseq = ''.join(new_subseq_codons)
                    if not 1 in [R in new_subseq for R in unamb_RE_seqs]:
                        seq_vec.append(new_subseq)
                        time = 0
                        for codon in new_subseq_codons:
                            time += speed_lookup[codon]
                        times_vec.append(time)
        best_option_index = times_vec.index(optimise_by[1](times_vec))
        new_subseq = seq_vec[best_option_index]
        new_seq = new_seq[:subseq_start] + new_subseq + new_seq[subseq_stop:]

    return new_seq

#==================================================================================================

def slow_down_seq(seq, codon_range = [], by=2, ref_table='Scer'):
    
    from CodOpY.analyse import time_seq,translate
    from CodOpY.misc import codon_choices
    
    #load auxiliary package data
    ref = load_from_Data(ref_table)
    time_dict_by_codon = dict(zip(ref['codon'], ref['decoding.time']))
    
    #add pseudo data for stop codons
    time_dict_by_codon['TGA'] = 0
    time_dict_by_codon['TAA'] = 0
    time_dict_by_codon['TAG'] = 0
    codon_seq = [seq[n:n+3] for n in range(0, len(seq),3)]
    if len(codon_range) == 2:
        codon_subseq = codon_seq[codon_range[0]:codon_range[1]+1]
    else:
        codon_subseq = codon_seq
    #determine existing and new target times for subseq
    is_time = time_seq(''.join(codon_subseq),ref_table=ref_table)['Average decoding time per codon']
    target_time = is_time * by
    #assemble a new sequence that is as close as possible to the new target time
    replace_codon_seq = []
    for codon in codon_subseq:
        aa = translate (codon)
        choice = codon_choices(aa)
        diffs = [abs(target_time - val) for val in list(choice['decoding.time'])]
        replace_codon_seq.append(choice.iloc[diffs.index(min(diffs))]['codon'])
    if len(codon_range) == 2:
        return_codon_seq = codon_seq[0:codon_range[0]] + replace_codon_seq + codon_seq[codon_range[1]+1:]
    else:
        return_codon_seq = replace_codon_seq
    return ''.join(return_codon_seq).replace('U','T')