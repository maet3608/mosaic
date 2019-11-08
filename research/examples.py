"""Runs some simple examples"""

import glob
import random
import pylab as pyl
from mods.seq import *
from mods.spectral import *
from mods.mosaic import *

if __name__ == "__main__":
    print "Started"

    filenames = glob.glob('data/MR_GR_new.ffa')
    #filenames = glob.glob('data/ABCA_new.ffa')
    #filenames = glob.glob('data/CCR2_CCR5_CCR3_new.ffa')
    #filenames = glob.glob('data/TLR1_TLR6_new.ffa')
    #filenames = glob.glob('data/Bactin_ARP1_ARP2_new.ffa')
    #filenames = glob.glob('data/*_new.ffa')
    for filename in filenames:
        print "loading ... "+filename
        #sequences = load_phylip(filename)
        sequences = load_fasta(filename)
        ngramize(sequences, 4)
        #sequences = randomize(sequences)

        #random.shuffle(sequences)                         # random sorting
        #sequences.sort(cmp=lambda a,b: len(a)-len(b))    # sort by length
        #sequences.sort(key=lambda s: s.name())           # sort by family
        sequences.sort(key=lambda s: s.letters())           # sort by sequence

        print "clustering ... "
        sequences,D,fi = sort_spectral_rec(sequences, 'L')     # sort recursively by spectrum

        print "analyzing..."
        draw_mosaic(filename,sequences,D,fi)
    pyl.show()
    print "finished"