import matplotlib

import glob
from re import match
from pylab import *
from scipy.stats import linregress
from random import random, shuffle
from mods.utils import *
from mods.seq import *



def pn(name):
    return name[:name.find('_')]
def ps(sequence):
    return pn(sequence.name())
  

def colordict(sequences):
    prefixes = enumerate(prefixset(sequences))
    return dict([(prefix,"bgrkcmy"[i%7]) for i,prefix in prefixes])

def prefixset(sequences):
    return set([ps(s) for s in sequences])


def read_D(filename):
    """reads an upper triangular matrix with phylogentic distances
       between species.
       Returns a dictionary with distances for species names."""
    D = {}
    f = open(filename)
    cols = f.readline().rstrip().split('\t')
    for line in f:
        elems = line.strip().split('\t')
        for i,e in enumerate(elems):
            if e!='' and i>0:
                D[(cols[i],elems[0])] = float(e)
                D[(elems[0],cols[i])] = float(e)
    f.close()
    return D


def dist(D,s1,s2):
    """Returns the phylogentic distance between two sequences"""
    return D[s1.name()[-2:],s2.name()[-2:]]


def optimizeD(D, sequences, name):
    ns = len(sequences) 
    c = []
    r = range(1,20)
    for n in r:
        for s in sequences: 
            s.calc_ngrams(n)
        x,y = [],[]
        for i in range(ns):
            for j in range(i+1, ns):
                s1 = sequences[i]
                s2 = sequences[j]
                x.append(dist(D,s1,s2))
                y.append(s1.distance(s2))
        corr = corrcoef(x,y)[0][1]                
        c.append(corr)
        #print n,corr
    #figure()    
    plot(r, c, '+-')
    #axis([1,20,0,1])
    xlabel("n-gram size")
    ylabel("correlation coefficient")
    #title(name)


def plotD(D, sequences, name, n):
    cd = colordict(sequences)
    ngramize(sequences,n)
    ns = len(sequences)  
    x,y,labels,colors = [],[],[],[]  
    for i in range(ns):
        for j in range(i+1, ns):
            s1 = sequences[i]
            s2 = sequences[j]
            x.append(dist(D,s1,s2))
            y.append(s1.distance(s2))
            labels.append(s1.name()[-2:]+s2.name()[-2:])
            colors.append(cd[pn(s1.name())])
    figure()    
    title(name)
    plot(x,y, 'o')
    #for i in range(len(labels)):
    #    text(x[i],y[i], labels[i], size=8)              
    #    text(x[i],y[i], labels[i], color=colors[i], size=8)
    (a,b,r,p,e) = linregress(x, y)
    print "a=%.2f b=%.2f, r=%.2f err=%.3f" % (a,b,r,e)

    xr = linspace(min(x), max(x))
    yr = polyval([a,b],xr)
    plot(xr,yr,'g-')        


def optimizeAUC(sequences, name):
    aucs = []
    r = range(1,20)
    for n in r:
        seqs1 = ngramize(sequences,n)
        seqs2 = ngramize(randomize(sequences),n)
        sims1 = similarities(seqs1, seqs1)
        sims2 = similarities(seqs1, seqs2)
        aucs.append(AUC(sims1+sims2,[1]*len(sims1)+[0]*len(sims2))) 
    #figure()    
    plot(r, aucs, '+-')
    axis([0, max(r), min(aucs)*0.9, max(aucs)*1.1])
    #title(name)
    xlabel("n")
    ylabel("AUC")


def optimizeFunc(sequences, name, f):
    m  = []
    r = range(1,20)
    for n in r:
        ngramize(sequences,n)
        sims = similarities(sequences, sequences)
        m.append(f(sims))
    #figure() 
    plot(r, m, '+-')
    xlabel('n')
    ylabel("std")
    #title(name)


def plotSimDistribution(sequences, filename, n):
    figure()
    ngramize(sequences,n)
    hist(similarities(sequences, sequences), 10)
    title(filename)


def std_entropy(sims, f):
    return std(map(lambda p : 0 if p<1e-8 else -p*log(p), sims)) 
    
 
def plotRegression():
    filename = 'data/MR_GR_new.ffa'
    sequences = load_fasta(filename) 
    D = read_D("data/D.mat")
    sequences = select_contains('AR', sequences)
    plotD(D, sequences, filename, 4)
        
    
def plotDivergenceCorrelation():
  filenames = glob.glob('data/MR_GR_new.ffa')
  #filenames = glob.glob('data/ABCA_new.ffa')
  #filenames = glob.glob('data/TLR1_TLR6_new.ffa')
  #filenames = glob.glob('data/Bactin_ARP1_ARP2_new.ffa')
  #filenames = glob.glob('data/*_new.ffa')
  labels = []
  for filename in filenames:
    print "loading ... "+filename
    sequences = load_fasta(filename)    
    D = read_D("data/D.mat")
    optimizeD(D, sequences, filename)
    for label in prefixset(sequences):
        labels.append(label)
        optimizeD(D, select_contains(label+'_', sequences), label)
  legend(labels, loc='best')  

  
def plotAUC():
  filenames = glob.glob('data/*_new.ffa')
  for filename in filenames:
    print "loading ... "+filename
    sequences = load_fasta(filename)    
    optimizeAUC(sequences, filename)
  labels = [match('data.(.+)_new',l).group(1) for l in filenames]  
  legend(labels, loc='best')  
 
      
def plotSTD():        
  filenames = glob.glob('data/*_new.ffa')
  for filename in filenames:
    print "loading ... "+filename
    sequences = load_fasta(filename)    
    optimizeFunc(sequences, filename, std)
  labels = [match('data.(.+)_new',l).group(1) for l in filenames]  
  legend(labels, loc='best')  

    
if __name__ == "__main__":
  print "Started"
  #plotRegression()
  #plotDivergenceCorrelation()
  #plotAUC()
  plotSTD()
  show()
  print "finished"
