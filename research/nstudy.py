import matplotlib

import glob
from re import match
from collections import defaultdict
from pylab import *
from random import shuffle, randint
from mods.utils import *
from mods.seq import *



def randomize(s):
    aas = [aa for aa in s.letters()]
    shuffle(aas)         # in-place
    return "".join(aas)


def artificial(length, as=20):
    return "".join([chr(64+randint(1,as)) for i in range(0,length)])


def ngram_freqs(letters, n):
    freqs = defaultdict(lambda: 0)
    num_ngrams = len(letters)-(n-1)
    for i in xrange(0, num_ngrams):
        freqs[letters[i:i+n]] +=1
    return freqs    

def avg_freq(letters, n):
    freqs = ngram_freqs(letters, n)
    return sum(f for f in freqs.values())/float(len(freqs.values()))
   
   
def max_freq(letters, n):
    freqs = ngram_freqs(letters, n)
    return max(f for f in freqs.values())


def avg_freq_seq(sequence, n):
    """Creates a set of n-grams"""
    #letters = sequence.letters()
    #letters = randomize(sequence)
    letters = artificial(3*len(sequence), 4)  #DNA
    #letters = artificial(len(sequence), 20)  #AA
    return avg_freq(letters,n)
    

          
      
def plotNumOverN():        
    filenames = glob.glob("data/*.ffa")
    freqsMat = []

    for filename in filenames:
        print "loading ... "+filename
        sequences = load_fasta(filename)  
    
        freqs = []
        n_sizes = range(1,10)
        for n in n_sizes:
            freq = sum(avg_freq_seq(s,n) for s in sequences)/float(len(sequences))
            freqs.append(freq)
        freqsMat.append(freqs)
        
        plot(n_sizes, freqs, "+-")
        #for i in range(0, len(n_sizes)): print("%3d %.3f"%(n_sizes[i],freqs[i]))
    
    freqsMat = array(freqsMat)
    print " ".join(["%5.3f"%float(val) for val in mean(freqsMat, 0)])
    xlabel("n-gram size")
    ylabel("average frequency")
    gca().set_xticks(n_sizes)
    gca().set_xticklabels(map(str, n_sizes))

    labels = [match('data.(.+)\.ffa',l).group(1) for l in filenames]  
    legend(labels, loc='best')  


def plotAvgFreqOverLen(): 
    letters = artificial(100000, 20)
    lens = range(10, 100000, 1000)
    freqs = [avg_freq(letters[0:l], 4) for l in lens]
    plot(lens, freqs, "-")    
    xlabel("sequence length")
    ylabel("average frequency")
    
    
    
    
if __name__ == "__main__":
  print "Started"
  plotNumOverN()
  #plotAvgFreqOverLen()
  show()
  print "finished"
