"""Draws a distance/similarity matrix as a network"""

import glob
from mods.seq import *
import networkx as nx
import pylab as pyl



def create_graph(name, sequences, cutoff=None):
  """Creates a graph based on the n-gram similarity of the sequences in the
  given list of sequences"""
  nodes = [s.name() for s in sequences]
  edges = [(s1.name(), s2.name(), s1.distance(s2)) for s1 in sequences for s2 in sequences]
  cutoff = cutoff if cutoff else pyl.median([d for (n1,n2,d) in edges])
  print "edge cutoff: %.2f"%cutoff
  edges = filter(lambda (n1,n2,d): d < cutoff, edges)
  G = nx.Graph(name=name)
  G.add_nodes_from(nodes)
  G.add_edges_from(edges)
  return G


def draw_labels(G, pos):
  colors = "bgrkcmy"
  prefixes = set([name[:-3] for name in G.nodes()])
  print "prefixes: "+",".join(prefixes)
  for i,prefix in enumerate(prefixes):
    draw_label(G, pos, prefix, colors[i%len(colors)])


def draw_label(G, pos, prefix, color):
  labels = dict([(n, str(n)) for n in G.nodes() if n.startswith(prefix)])
  nx.draw_networkx_labels(G, pos, labels=labels, font_color=color, font_size=8)


def draw_edges(G, pos, cutoff=None):
  if cutoff:
    print "edge drawing cutoff: %.2f"%cutoff
    edgelist = filter(lambda (n1,n2,d): d < cutoff, G.edges_iter(data=True))
    nx.draw_networkx_edges(G, pos, edgelist=edgelist, edge_color="b", alpha=0.2)
  else:  
    edge_color = [d for (n1,n2,d) in G.edges_iter(data=True)]
    nx.draw_networkx_edges(G, pos, edge_color=edge_color, edge_cmap = pyl.cm.gray_r)


def draw_nodes(G, pos):
  nx.draw_networkx_nodes(G, pos, node_size = 1000, node_color="white")
  #nx.draw_networkx_nodes(G,pos)


def draw_graph(G):
  pyl.figure()
  #pos = nx.circular_layout(G)
  pos = nx.spring_layout(G)
  #pos = nx.spectral_layout(G)
  draw_edges(G, pos, None)
  draw_nodes(G, pos)
  draw_labels(G, pos)
  pyl.title(G)
  #pyl.axis('off')
  #pyl.savefig("graph.png")
  


if __name__ == "__main__":
  print "Started"
  for filename in glob.glob('data/MR_GR_new.ffa'):
  #for filename in glob.glob('data/*.ffa'):
    print "loading ... "+filename
    sequences = load_fasta(filename)
    ngramize(sequences, 4)
    print "analyzing..."
    G = create_graph(filename, sequences, None)
    draw_graph(G)
  pyl.show()
  print "finished"
