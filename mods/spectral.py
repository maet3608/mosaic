"""Spectral rearrangment."""

try: import Psycho
except ImportError: pass 

import glob
from pylab import *
from utils import *
from seq import *


def laplacian(A, Ltype='L'):
    """Calculates a Laplacian matrix.
    A - Affinity matrix.
    Ltype - type of Laplacian: L, Lnorm, Lsym
    """
    A[A<0.01] = 0                 # Remove small similarities => sparse A
    I = identity(len(A))          # Identity matrix
    s = A.sum(axis=0)             # row sum

    if Ltype=='L':                # unnormalized Laplacian
        D = diag(s)               # D
        L = D-A                   # D-A
    elif Ltype=='Lnorm':          # normalized Laplacian
        D = diag(1/s)             # D^-1 
        L = I - dot(D,A)          # L = I - D^-1 A
    elif Ltype=='Lsym':           # normalized symmetric Laplacian
        D = diag(s**(-0.5))       # D^-1/2 
        L = I - dot(dot(D,A),D)   # L = I - D^-1/2 A D^-1/2 : Shi-Mali
    else:
        raise "Invalid matrix type: "+Ltype
      
    return L


def fiedler(A, Ltype='L'):
    """Calculates the fiedler vector
    A - Affinity matrix.
    Ltype - type of Laplacian: L, Lnorm, Lsym
    """
    L = laplacian(A, Ltype)
    w,v = eigh(L)               # Eigenvalues/vectors of the laplacian 
    si = argsort(w)             # sorting index acc. to eigenvalues

    ev = v[:,si[1]]             # fiedler vector
    fi = argsort(ev)            # sorting index acc. to fiedler vector

    #f = ev[fi]
    #figure(); plot(ev[fi], "o-b");title('Fiedler vector')
    #axis([-1, len(f), min(f)*1.1, max(f)*1.1])
    #figure(); plot(w[si], "x-r"); title('eigenvalues')
        
    #figure()
    #matshow(v, cmap=cm.Greens); title('Eigenvectors')
    #matshow(v[:,si][fi,:], cmap=cm.Greens); title('Eigenvectors sorted')
    #matshow(A, cmap=cm.Blues); title('Affinity matrix')
    #matshow(A[:,fi][fi,:], cmap=cm.Blues); title('Affinity matrix sorted')
    return ev


def affinity_matrix(sequences, sigma=0):
    """Returns the affinity matrix and the distance matrix
    sequences - list of sequences D is computed from
    sigma - radius of Gaussian kernel for affinity matrix, 0 means automatic 
    """
    n = len(sequences)
    D = zeros((n,n))
    for i in xrange(n):
        for j in xrange(i+1,n):
            D[i][j] = D[j][i] = sequences[i].distance(sequences[j])           
    sigma2 =  sigma if sigma>0 else 2*mean(D)**2
    return exp(-D**2/sigma2),D


def min_cut(ev):
    """Returns the index where the elements of the Fiedler vector change
    it sign and the difference between the two elements
    ev -- fiedler vector 
    """
    for i in range(1,len(ev)):
        if ev[i-1] < 0 < ev[i]: return (i, ev[i]-ev[i-1])
    return (None,0)   # can't cut, one uniform cluster


def outliers(ev, outlier=3.0):
    """Returns a list of indices that indicates the elements of the
    fiedler vector that are NOT outliers, with an outlier being defined
    as an element more then outlier standard deviations away from the mean
    ev -- fiedler vector 
    outlier -- number of standard deviations that define an outlier
    """
    t = outlier*std(ev)
    v = abs(ev-mean(ev))
    idx = arange(len(ev))
    print "outliers:",idx[v>t]
    return idx[v<t]


def filter_outliers(A, outlier=3.0, Ltype='L'):
    """Filters all outliers from the affinity matrix A by analyzing the
    standard deviation of the fiedler vector
    A -- affinity matrix
    outlier -- threshold for outliers
    Ltype -- type of Laplacian
    """
    ev = fiedler(A, Ltype)
    oi = outliers(ev, outlier)
    return oi, A[:,oi][oi,:]


def rec_reordering(A, o, Ltype='L'):
    """Recursive spectral reordering.
    A - affinity matrix
    o - row/col ordering of A
    Ltype - type of Laplacian: L, Lnorm, Lsym
    Returns a reordering index for A
    """
    if len(o) < 2: return o         # nothing left to reorder
    L = laplacian(A, Ltype)         # calc laplacian
    w,v = eigh(L)                   # Eigenvalues/vectors of the laplacian
    ev = v[:,argsort(w)[1]]         # fiedler vector
    fi = argsort(ev)                # sorting index acc. to fiedler vector

    if list(fi)==range(len(fi)):    # no change in ordering
        return o                    # finish recursion
    
    A = A[:,fi][fi,:]               # reorder matrix
    o = o[fi]                       # reorder ordering
    
    # split matrix at min cut and call recursivly for sub matrices
    mc,d = min_cut(ev[fi])          # find min cut
    if not mc: return o             # can't split
    
    o1 = rec_reordering(A[:mc,:mc], o[:mc], Ltype)
    o2 = rec_reordering(A[mc:,mc:], o[mc:], Ltype)
    return concatenate((o1,o2))


def sort_spectral(sequences, Ltype='L', sigma=0):
    """Returns a list of sequences sorted according to the components of
    the fiedler vector resulting from a spectral rearrangment
    sequences - list of sequences to sort
    Ltype - type of laplacian
    sigma - RBF radius for affinity matrix
    """
    A,D = affinity_matrix(sequences, sigma)
    ev = fiedler(A, Ltype)
    return [s for s,f in sorted(zip(sequences, ev), key=lambda t:t[1])]


def sort_spectral_rec(sequences, Ltype='L', sigma=0):
    """Returns a list of sequences sorted according to the components of
    the fiedler vector resulting from a *recursive* spectral rearrangment
    sequences - list of sequences to sort
    Ltype - type of laplacian
    sigma - RBF radius for affinity matrix
    """
    A,D = affinity_matrix(sequences, sigma)
         
    o = rec_reordering(A, array(range(len(A))), Ltype)
    fi = fiedler(A[:,o][o,:], Ltype)
    sequences = [sequences[i] for i in o]
        
    return (sequences, D[:,o][o,:], fi)


def load_matrix(filepath):
    """Loads a (distance) matrix.
    filepath -- path to matrix to load
    """
    f = open(filepath)
    mat = []
    for line in f:
        mat.append([float(e) for e in line.strip().split('\t')[1:]])
    f.close()
    return array(mat)



if __name__ == "__main__":
    """Just a usage example"""
    filenames = glob.glob('data/MR_GR.ffa')
    for filename in filenames:
        print "loading ... "+filename
        sequences = load_fasta(filename)
        ngramize(sequences, 4)
        sequences.sort(key=lambda s: s.name())
        sort_spectral_rec(sequences)
        #sort_spectral(sequences)
    show()
