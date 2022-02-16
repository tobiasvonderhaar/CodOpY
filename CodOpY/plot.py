import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
from CodOpY.optimise import opt_seq, translate
from CodOpY.misc import load_from_Data

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
