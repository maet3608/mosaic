"""Calculation and visualization of n-gram profile"""

try: import Psycho
except ImportError: pass 

from pylab import *
from seq import *
from utils import *
from collections import defaultdict



def calc_profile(sequences, n=4):
    """Calculates an n-gram profile over sequences,
    which is the frequency distribution of n-grams and
    an averaged position.
    sequences -- list of sequences
    n -- n-gram size
    """
    profile = defaultdict(lambda: [0,0]) 
    for sequence in sequences:
        letters = sequence.letters()
        for i in xrange(len(letters)-(n-1)):
            ngram = letters[i:i+n]
            info = profile[ngram]
            info[0] += 1                   # frequency
            info[1] += i                   # positions
    for ngram, info in profile.items():
        info[1] = info[1]/info[0]   #average position
    return profile


def calc_mprofile(profile):
    """Creates a profile that contains only the maximum frequencies over
    positions
    profile -- a profile created by calc_profile()
    """
    mprofile = defaultdict(lambda: 0)
    for f,p in profile.values():
        mprofile[p] = max(mprofile[p],f)
    return mprofile        


def print_profile(profile, show=10):
    """Prints out the profile n-grams sorted according to position.
    profile -- profile to print
    show -- number of n-grams to print""" 
    items = sorted(profile.items(), key=lambda t: -t[1][0])
    print("Number of n-grams: %d"%len(items))
    for ngram, info in items[:show]:
        print("%s  %d %d"%(ngram, info[0],info[1]))


def plot_mprofile(sequences, name, n=4):
    """Plots a maximum profile.
    sequences -- list of sequences
    name -- plot title
    n -- size of n-grams
    """
    mprofile = calc_mprofile(calc_profile(sequences, n))
    items = sorted(mprofile.items(), key=lambda t: t[0])
    p,f = zip(*items)  # positions and frequencies
    f = mean_filter(f,len(f)/10)
    plot(p,f,'r-')
    title(name)


def plot_profile(profile, name):
    """Plots a  profile.
    profile - profile created by calc_profile()
    name -- plot title
    """
    values = profile.values()
    f,p = zip(*values)  # frequencies and positions
    plot(p,f,'.')
    title(name)

 


def mean_filter(x, k):
    """As the name say, a mean filter.
    x - array to filter
    k - window size.
    """
    l = len(x)
    y = [0]*l
    xs = [x[0]]*ceil(k/2)
    xe = [x[l-1]]*floor(k/2)
    x = xs + list(x) + xe
    for i in range(l):
        y[i] = mean(x[i:i+k])
    return y        

  
  
if __name__ == "__main__":
    """Just a usage example"""
    
    filenames = glob.glob('data/MR_GR.ffa')
    #filenames = glob.glob('data/ABCA.ffa')
    #filenames = glob.glob('data/TLR1_TLR6.ffa')
    #filenames = glob.glob('data/Bactin_ARP1_ARP2.ffa')
    #filenames = glob.glob('data/*_new.ffa')
    for filename in filenames:
        sequences = load_sequences(filename)
        #profile = calc_profile(select_contains('MR_', sequences))
        profile = calc_profile(sequences)
        figure(); plot_profile(profile, filename)
        figure(); plot_mprofile(calc_mprofile(profile), filename)
        #print_profile(profile)
    show()
