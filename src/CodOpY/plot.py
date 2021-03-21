def count_aas(seq):
    '''Produces a bar graph of amino acid counts in a protein sequence'''

    import matplotlib.pyplot as plt

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
