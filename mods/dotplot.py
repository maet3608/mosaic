"""Shows a dot plot based on n-gram matches between two sequences."""

try: import Psycho
except ImportError: pass 

from seq import *
from pylab import *
from utils import *



def calc_ps(s1, s2, n):
    """Calculates a list of n-gram matches between the two sequences.
    s1 -- first sequence
    s2 -- second sequence
    n -- n-gram size
    """
    n1 = len(s1)-n-1
    ps = []
    for i in xrange(n1):
        ngram = s1[i:i+n]
        j = s2.find(ngram)  
        while j >= 0:
            ps.append((i,j))
            j = s2.find(ngram, j+1)
    return ps


def mat_create(s1, s2, n, c):
    """Creates the dot matrix. Each dot indicates a matching n-mer.
    s1 - first sequence
    s2 - second sequence
    n - size of n-gram, e.g. 4
    c - scaling factor. Summarizes dots within cells of the given size.
    """
    n1, n2 = len(s1)-n-1, len(s2)-n-1
    mat = zeros((ceil(n1/float(c))+1, ceil(n2/float(c))+1))
    ps = []
    for i in xrange(n1):
        ngram = s1[i:i+n]
        j = s2.find(ngram)  
        while j >= 0:
            mat[i/c,j/c] += 1
            ps.append((i,j))
            j = s2.find(ngram, j+1)
    return mat


def plot_dot(s1, s2, n=4, c=20, bin=False):
    """Creates a dot plot.
    s1 - first sequence
    s2 - second sequence
    n - size of n-gram, e.g. 4
    bin - use only two colors
    """
    s1.calc_ngrams(n)
    s2.calc_ngrams(n)
    s = s1.similarity(s2)
    mat = mat_create(s1,s2,n, c)
    if bin : mat[mat>0] = 1
    nr, nc = mat.shape
    extent = [-0.5, nc-0.5, nr-0.5, -0.5]
    imshow(mat, extent=extent, cmap=cm.Blues, origin='upper', interpolation='nearest')
    text(0.05*len(s2)/c,0.95*len(s1)/c, "s=%.2f, n=%d, c=%d"%(s,n,c), color='blue') 
    xlabel(s1.name())   
    ylabel(s2.name())   
    title("Dot plot:  "+s1.name()+':'+s2.name())


def plot_dots(sequences, n=4, c=20):
    """Creates dot plots for all the given sequences."""
    ns = len(sequences)
    for i in xrange(ns):
        for j in xrange(i,ns):            
            ax = subplot(ns,ns, i*ns+j+1)
            ax.axis('off')  
            plot_dot(sequences[i],sequences[j],n,c)


def plot_dot_multi(sequences, n=4, c=5):
    """Creates dot plots for all the given sequences."""
    def fill_mat(mat, s1, s2, n, c):
        for i in xrange(len(s1)):
            ngram = s1[i:i+n]
            j = s2.find(ngram)  
            while j >= 0:
                mat[i/c,j/c] += 1
                j = s2.find(ngram, j+1)
    ml = max(map(len, sequences)) #longest sequence
    rc = ceil(ml/float(c))
    mat = zeros((rc,rc))
    ns = len(sequences)
    for i in xrange(ns):
        for j in xrange(i+1,ns):            
            fill_mat(mat, sequences[i], sequences[j], n, c)
    extent = [-0.5, ml-0.5, ml-0.5, -0.5]
    imshow(mat, extent=extent, cmap=cm.hot, origin='upper', interpolation='nearest')




  
if __name__ == "__main__":
    """Just a usage example"""
    
    sequences = load_fasta("data/MR_GR.ffa")
    #c = max(map(len, sequences))/30
   
    plot_dot(take('AR_Hs',sequences),take('AR_Hs',sequences),n=4,c=10,bin=True)
    #plot_dot(sequences[0],sequences[200],n=2,c=2,bin=True)
    #plot_dots([sequences[i] for i in map(int, linspace(0, len(sequences), num=6, endpoint=False))], c=c)
    #plot_dots(sequences, c=c)
    #plot_dot_multi(select_starts("MR_",sequences)+ select_starts("GR_",sequences), c=5)
    #map(mod_alpha, sequences)
    #plot_dot_multi(sequences, n=3, c=5)
    #plot_dot_multi(sequences)
                   
    show()

