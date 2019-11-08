"""The main application that reads data files and plots MOSAIC matrix"""

try: import Psycho
except ImportError: pass 

import pylab as pyl
from mods.dotplot import plot_dot
from mods.profile import plot_mprofile
from mods.utils import *
from mods.seq import *
from mods.spectral import *



def show_flowers():
    """Just some flowers"""
    im = matplotlib.image.imread("doc/pic/flowers.png")
    fig = figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.imshow(im, interpolation='bilinear')


def draw_mosaic(name, sequences, D, fi, options): 
    """draws the MOSAIC plot.
    name -- name of data set
    sequences -- list of sequences
    D -- distance matrix
    options -- command line options
    """  
    class Var: pass      #just because handling of global var is ugly otherwise
    Var.doPlot = False   #dot and profile plots only when shift key pressed
    si,sj = None,None    #start coord. of selection rectangle
    
    def indices(event):
        if event.xdata == None or event.ydata == None: return (None,None)
        return int(event.xdata+0.5), int(event.ydata+0.5)
    def onKeyPressed(event):
        Var.doPlot = event.key == "control"
        fig.canvas.mpl_disconnect(kpi)
    def onKeyReleased(event):
        Var.doPlot = False
        global kpi
        kpi = fig.canvas.mpl_connect('key_press_event', onKeyPressed)
    def onButtonPressed(event):
        global si,sj
        if not Var.doPlot: return
        Var.doPlot = False
        si,sj = indices(event)
    def onButtonReleased(event):
        if not Var.doPlot: return
        Var.doPlot = False
        ei,ej = indices(event)
        if (ei,ej) == (None,None): return
        if (ei,ej) == (si,sj):
            figure()
            plot_dot(sequences[ei], sequences[ej], n=options.n, 
                     c=options.c, bin=options.binary);
            show()
        else:
            i1, i2 = min([ei,ej,si,sj]), max([ei,ej,si,sj])
            figure()
            plot_mprofile(sequences[i1:i2+1], "%d-%d"%(i1,i2), options.n)
            show()
            
    if not options.no_mosaic:
        n = len(sequences)
        cd = colordict(sequences)
        fig = figure(figsize=(8.0,8.2))
        matshow(D, fignum=fig.number, cmap=cm.hot, aspect='equal')
        title("Mosaic: "+name)
        fs = options.fontsize
        for i,label in enumerate(s.name() for s in sequences):
            color = "k" if options.black_labels else cd[pn(label)]
            text(i, n+0.1, label, size=fs, color=color, rotation='vertical', ha='center', va='top')
            text(n+0.1, i, label, size=fs, color=color, rotation='horizontal', ha='left', va='center')
        fig.canvas.mpl_connect('button_release_event', onButtonReleased)
        fig.canvas.mpl_connect('button_press_event', onButtonPressed)
        fig.canvas.mpl_connect('key_release_event', onKeyReleased)
        kpi = fig.canvas.mpl_connect('key_press_event', onKeyPressed)
   
    if not options.no_fiedler:
        figure()
        plot(fi, "o-b")
        axis([-1, len(fi), min(fi)*1.1, max(fi)*1.1])
        title("Fiedler: "+name)
        
    if options.flowers:
        show_flowers()
     
     
def create_parser():
    """Creates the command line parser"""
    from optparse import OptionParser, SUPPRESS_HELP 

    parser = OptionParser(usage="usage: %prog [options] filepath", 
                          version="mosaic 1.01 (15.10.2009)")
                          
    parser.add_option("-f", "--format", dest="format", default="FASTA",
                        help="input format: FASTA or PHYLIP")
    parser.add_option("-n", "--ngram_size", dest="n", type="int", default=4,
                        help="n-gram size")
    parser.add_option("-t", "--type", dest="type", default="L",
                        help="type of the Laplacian: L, Lnorm, Lsym")
    parser.add_option("-r", "--radius", dest="sigma", type="float", default=0.0,
                        help="radius of Gaussian kernel")
    parser.add_option("-s", "--font_size", dest="fontsize", type="int", default=9,
                        help="font size for labels of mosaic plot")
    parser.add_option("-k", "--black_labels",   
                        action="store_true", dest="black_labels", default=False,
                        help="black labels of mosaic plot (otherwise automatic coloring)")
    parser.add_option("-c", "--cell_size", dest="c", type="int", default=15,
                        help="cell/grid size of dot plot")
    parser.add_option("-b", "--binary",
                        action="store_true", dest="binary", default=False,
                        help="switches to black and white colors for the dot plot")
    parser.add_option("-i", "--no_fiedler", 
                        action="store_true", dest="no_fiedler", default=False,
                        help="suppress fiedler vector")
    parser.add_option("-m", "--no_mosaic", 
                        action="store_true", dest="no_mosaic", default=False,
                        help="suppress mosaic plot")
    parser.add_option("-v", "--verbose",
                        action="store_true", dest="verbose", default=False,
                        help="print status messages to stdout")
    parser.add_option("-w", "--flowers",
                        action="store_true", dest="flowers", default=False,
                        help=SUPPRESS_HELP)     
    return parser
    
    
    
if __name__ == "__main__":
    """The main method"""
    from glob import glob
     
    parser = create_parser()
    (options, args) = parser.parse_args()
    
    if len(args) != 1: 
        parser.error("Filepath is missing! Use -h for help.")

    filepath = args[0]
    verbose, isFasta = options.verbose, options.format == "FASTA"
    
    if verbose: print "running ... "
    filenames = glob(filepath)
    if not filenames:
        parser.error("Can't find files under: "+filepath)
        
    for filename in filenames:
        if verbose: print "loading ... "+filename
        sequences = load_fasta(filename) if isFasta else load_phylip(filename)
        ngramize(sequences, options.n)
        sequences.sort(key=lambda s: len(s))          
 
        if verbose: print "analyzing ... "
        sequences,D,fi = sort_spectral_rec(sequences, options.type, options.sigma)    

        if verbose: print "drawing ..."
        draw_mosaic(filename, sequences, D, fi, options)
    pyl.show()   
     
    if verbose: print "finished"
    
  
