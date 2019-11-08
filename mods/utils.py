import matplotlib

import glob
from pylab import median, mean, std
from random import shuffle
from seq import *


def pn(name):
    """Returns the prefix of a name, with the prefix
    separated by '_' for the remainder of the name.
    name - string
    """
    return name[:name.find('_')]
    
    
def ps(sequence):
    """Returns the prefix of a sequence name, with the prefix
    separated by '_' for the remainder of the name.
    sequence - sequence the name is extracted from.
    """
    return pn(sequence.name())
  

def colordict(sequences):
    """ returns a dictionary that maps prefixes of sequence names to
    label colors for the dot plot"""
    prefixes = enumerate(set([ps(s) for s in sequences]))
    return dict([(prefix,"bgrkcmy"[i%7]) for i,prefix in prefixes])


def mod_alpha(sequence):
    """Modifies the alphabet and replaces the letters of a sequence
    by letters drawn from a different alphabet"""
    code = { 'H':'0','R':'0','K':'0',
            'D':'1','E':'1','N':'1','Q':'1',
            'S':'2','T':'2','P':'2','A':'2','G':'2',
            'M':'3','I':'3','L':'3','V':'3',
            'F':'4','Y':'4','W':'4',
            'C':'5',
            'X':'6' }
    sequence._letters = "".join(code[aa] for aa in sequence._letters)     


def similarities(seqs1, seqs2):
    """Calculates the similarities between two lists of sequences.
    seqs1 -- first list of sequences
    seqs2 -- second list of sequences
    """
    sims = []
    ns = len(seqs1)
    for i in range(ns):
        for j in range(i+1, ns):
            s1 = seqs1[i]
            s2 = seqs2[j]
            sims.append(s1.similarity(s2))
    return sims

  
def select_starts(prefix, sequences):
    """Selects sequences from the given list of sequences which info attribute
    start with the given prefix."""
    print("select sequences: "+prefix)
    return [s for s in sequences if s.info().startswith(prefix)]
    
    
def select_contains(substr, sequences):
    """Selects sequences from the given list of sequences which name
    contain the given substring."""
    print("get sequences: "+substr)
    return [s for s in sequences if substr in s.name()]
    
    
def take(name, sequences):
    """Takes a sequence with the given name from the list of sequences"""
    for s in sequences:
        if s.name() == name: return s
    return None    


def randomize(sequences):
    """Glues all sequences together, randomly shuffles the letters
       and the splits the shuffled letter sequence in sub-sequences
       of the same lengths as the original sequences."""
    print "randomizing sequences"
    aas = []
    for s in sequences:
        aas += [aa for aa in s.letters()]
    shuffle(aas)
    rnd_sequences = []
    i = 0
    for s in sequences:
        n = len(s)
        letters = "".join(aas[i:i+n])
        i += n
        rnd_sequences.append(Seq('R'+s.name(), letters))
    return rnd_sequences


def threshold(sequences, n):
    """Calculates a cutoff threshold for the similarity between
       sequences"""
    sequences = ngramize(randomize(sequences),n)
    return(max(similarities(sequences, sequences)))
    
    
def AUC(output, target, classID1=0, classID2=1):
    """Area Under the ROC curve. 0.5 means no correlation.
       1.0 means perfect correlation.
       output - classifier output. Real value that represents confidence.
       target - target output.
       classID1 - ID of the first class
       classID2 - ID of the second class.
    """
    idx = sort_index(output)     # sorting index according to output
    np = target.count(classID1)  # number positives
    nn = target.count(classID2)  # number negatives
    tp, fp = 0.0, 0.0            # true and false positives
    ox, oy = 0.0, 0.0            # old pos of ROC curve
    old_i = None                 # old index
    auc = 0.0                    # AUC
    for i in idx:
        if target[i] == classID1: tp += 1.0
        elif target[i] == classID2: fp += 1.0
        else : continue

        if old_i and not output[old_i] == output[i]:
            x,y = fp/nn, tp/np
            auc += (x-ox)*(y+oy)/2.0
            ox, oy = x,y
        old_i = i
    auc += (1-ox)*(1+oy)/2.0
    return auc if auc > 0.5 else 1.0-auc

    
def sort_index(v):
    """Returns a sorting index for the elements in v.
       The smallest element in v has index zero
    """
    return [i for i,e in sorted(enumerate(v), key=lambda t: t[1])]

    
    
if __name__ == "__main__":
    """Just a usage example"""
    for filename in glob.glob('data/MR_GR.ffa'):
        print "loading ... "+filename
        sequences = load_sequences(filename, 0)
        sequences = select(sequences, 'MR_')
        sequences =  randomize(sequences)
        for s in sequences:
            print s
