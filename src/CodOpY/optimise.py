def retrieve_kazusa(taxID):
    '''Returns a table with codon usage frequencies (per 1000 nt) from the Kazusa website'''

    from CodOpY.misc import initiate_code_table
    import requests
    import pandas as pd

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

def opt_seq(seq,diversify=[],diversify_range=0.2,ref_table = 'Scer',optimise_by=['decoding.time',min]):

    import pandas as pd
    import random



    #determine which ref_table to use. If this is the default string 'Scer',
    #load this from the package Data
    if ref_table == 'Scer':
        #prepare package data for use
        try:
            import importlib.resources as pkg_resources
        except ImportError:
            # Try backported to PY<37 `importlib_resources`.
            import importlib_resources as pkg_resources
        from . import Data  # relative-import the *package* containing the data
        #import the stored data for S cerevisiae
        with pkg_resources.open_text(Data, 'Scer.csv') as read_file:
            parameterset = pd.read_csv(read_file)

    #define the reverse translation dictionary: which codons should be considered for which amino acid?
    #if an amino acid is not in the diversify list, use the fastest codon
    #if an amino acid is in the diversify list, diversify codon choice by random draw from all codons for which the diversify_range threshold applies
    reverse_dict = {}
    diversify_threshold = optimise_by[1]([max(codons[optimise_by[0]]*diversify_range),max(codons[optimise_by[0]]*(1- diversify_range))])
    for aa in codons['one.letter'].unique():
        codons_for_aa = codons.loc[codons['one.letter']==aa]
        if aa in diversify:
            acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']<diversify_threshold]['codon'])
        else:
            min_decoding_time_for_aa = min(codons_for_aa['decoding.time'])
            acceptable_codons_for_aa = list(codons_for_aa.loc[codons_for_aa['decoding.time']==min_decoding_time_for_aa]['codon'])
        reverse_dict[aa] = acceptable_codons_for_aa
    codon_seq = []
    for seq_aa in seq:
        codon_seq.append(random.choice(reverse_dict[seq_aa]))
    return ''.join(codon_seq).replace('U','T')


#==================================================================================================
