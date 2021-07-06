import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import requests
from CodOpY.misc import initiate_code_table

def retrieve_kazusa(taxID):
    '''Returns a table with codon usage frequencies per 1000 nt from the
    Kazusa website.

    Parameters
    ==========
    taxID : int or str
        taxID can be the NCBI Taxonomy ID of an organism or the latin name (or
        part thereof) of an organism.

    Returns
    =======
    pandas.core.frame.DataFrame
        A dataframe containing codons, amino acid abbreviations and the
        relative usage frequency per 1000 codons for each of the 64 possible
        RNA codons.
    '''

    if type(taxID) == str:
        if not taxID.isdecimal():
            search_result = requests.get('http://www.kazusa.or.jp/codon/cgi-bin/spsearch.cgi?species=' + taxID.replace(' ','+') + '&c=i')
            if '\nNot found\n' in search_result.text:
                return 'Search term not found'
            result_lines = search_result.text.split('<A')
            ids = []
            print('Available entries include:\n')
            for line in result_lines[:-1]:
                if line[1:5] == 'HREF':
                    idx = line.split('=')[2].split('\"')[0]
                    species = line.split('<I>')[1].split('</I>')[0]
                    if idx.isdecimal():
                        ids.append(idx)
                        print(species + ': ' + idx)
            if len(ids) == 1:
                print('Using unique result')
                taxID = ids[0]
            else:
                query = input('\nPlease select a numerical ID')
                if query not in ids:
                    return 'Invalid ID'
                else:
                    taxID = query

    base_frame = initiate_code_table()
    base_frame['codon'] = base_frame['codon'].str.replace('T','U')
    search_result = requests.get('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=' + str(taxID))
    if 'Not found' in search_result.text:
                return 'Search term not found'
    species = search_result.text.split('\n')[7].split('<i>')[1].split(' </i>')[0]
    print('\nRetrieving data for ' + species)
    if len(search_result.text.split('PRE')) == 1:
        print('no codon usage data found for this taxonomy ID')
        return
    else:
        result_table = search_result.text.split('PRE')[1].replace('\n','').replace('>','').replace('</','')
        result_table = result_table.split(')')
    #retrieve information from lines
    codons,frequency = [],[]
    for line in result_table:
        #remove leading spaces
        if len(line)>1:
            while line[0] == ' ':
                line=line[1:]
            codons.append(line[:3])
            frequency.append(float(line[line.find(' ')+1:line.find('(')]))
    results_frame = pd.DataFrame({'codon':codons,'usage.frequency':frequency})
    results_frame = base_frame.merge(results_frame, how = 'outer',on='codon')
    return results_frame.sort_values(by='codon').reset_index(drop=True)

#==================================================================================================

def opt_seq(seq,diversify=['K','N','I','H','V','G','D','Y','C','F'],
            diversify_range=0.2,ref_table = 'Scer',
            optimise_by=['decoding.time',min]):

    '''
    Makes an optimised DNA sequence corresponding to an input amino acid sequence.

    Parameters
    ==========
    seq : str
        The amino aid sequence to be otimised

    diversify : list of str
        A list specifying individual amino acids for which codons should be
        diversified.

    diversify_range : float
        The proportion over which the optimisation parameters that is allowed to
        vary for diversified amno acids, relative to 1.

    ref_table : str
        the name of the data file containing the codon data.

    optimise_by : list of str and function
        The str part of optimise_by specifies which cloumn of the ref_tble shold be
        used for optimisation. Function can be min or max and specifies whether the
        highestt or lowest value should be selected for optimisation.

    Returns
    =======
    str
        Returns a DNA sequence string.
    '''

    #prepare package data for use
    try:
        import importlib.resources as pkg_resources
    except ImportError:
        # Try backported to PY<37 `importlib_resources`.
        import importlib_resources as pkg_resources
    from . import Data  # relative-import the *package* containing the data
    #import the stored data for the dataset in question
    try:
        with open(Data.__path__[0] + '/' + ref_table + '.csv') as read_file:
            parameterset = pd.read_csv(read_file)
    except:
        raise ValueError('ref_table does not refer to a valid dataset.')

    #define the reverse translation dictionary: which codons should be considered for which amino acid?
    #if an amino acid is not in the diversify list, use the fastest codon
    #if an amino acid is in the diversify list, diversify codon choice by random draw from all codons for which the diversify_range threshold applies
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
                return
        else:
            best_decoding_time_for_aa = optimise_by[1](codons_for_aa['decoding.time'])
            acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']==best_decoding_time_for_aa]['codon'])
        reverse_dict[aa] = acceptable_codons_for_aa
    codon_seq = []
    for seq_aa in seq:
        codon_seq.append(random.choice(reverse_dict[seq_aa]))
    return ''.join(codon_seq).replace('U','T')


#==================================================================================================

def remove_RE(site, test_seq, ref_table = 'Scer',optimise_by=['decoding.time',min],suppress_not_found=False):

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
    with open(Data.__path__[0] + '/' + ref_table + '.csv') as read_file:
        codons = pd.read_csv(read_file)

     #Scer = pd.read_csv('src/CodOpY/Data/' + ref_table + '.csv')

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
#================================================================================

def test_RE(REs, test_seq):
    '''
    Tests whether REs are present in a sequence.

    Parameters
    ==========
    REs : str or list of str
        The name, or list of names, of the enzymes to be tested.

    test_seq : str
        The DNA sequence to be tested.

    Returns
    =======
    list
        A list of names of those enzymes for which sites were found in
        test_seq, or an empty list of none of the enzyme sites was found.
    '''
    if type(REs) != list:
        REs = [REs]

    found_REs = []
    for RE in REs:
        test = remove_RE(RE, test_seq, suppress_not_found=True)
        if type(test) != str:
            pass
        else:
            found_REs.append(RE)
    return found_REs

#================================================================================

def translate(seq):
    '''
    Returns an amino acid sequenc etranslated for na input DNA sequence.

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
    try:
        import importlib.resources as pkg_resources
    except ImportError:
        # Try backported to PY<37 `importlib_resources`.
        import importlib_resources as pkg_resources
    from . import Data  # relative-import the *package* containing the data
    #import the stored data for S cerevisiae
    with open(Data.__path__[0] + '\\' +  ref_table + '.csv') as read_file:
        ref = pd.read_csv(read_file)
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
