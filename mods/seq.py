"""Describes biological sequences, e.g. DNA"""

try: import Psycho
except ImportError: pass 



class Seq(object):
  """Represents a biological sequence"""
  
  def __init__(self, name, letters):
    """Constructor:
    name -- Name of the sequence
    letters -- String with the sequence letters"""
    self._name    = name
    self._letters = letters.upper()
    self._ngrams  = None

  def __getitem__(self, key):
    return self._letters[key]
  
  def calc_ngrams(self, n):
    self._ngrams = ngrams(self._letters, n,n)

  def ngrams(self):
    """Getter for the set of ngrams."""
    return self._ngrams

  def similarity(self, seq):
    # max -> global alignment, min -> local alignment
    c = float(min(len(seq.ngrams()), len(self._ngrams)))
    return len(seq.ngrams().intersection(self._ngrams))/c if c else 0

  def distance(self, seq):
    return 1.0-self.similarity(seq)

  def name(self):
    return self._name

  def letters(self):
    """Getter for the sequence letters"""
    return self._letters

  def find(self, ngram, start=0):
    """Finds the given ngram within the sequence. Returns the
    index of the match or -1 if no match could be found."""
    return self._letters.find(ngram, start)

  def __len__(self):
    return len(self._letters)
  
  def __str__(self):
    """String representation of the sequence"""
    return "%s: %s"%(self._name, self._letters)


def ngrams(letters, n1, n2):
  """Creates a set of n-grams ranging from n1 to n2
  letters -- string of letters
  n1 -- minimum n-gram size
  n2 -- maximm n-gram size 
  """
  ngrams = set()
  for n in xrange(n1, n2+1):
    for i in xrange(0, len(letters)-(n-1)):
      ngrams.add(letters[i:i+n])
  return ngrams
  

def ngramize(sequences, n):
  """Creates the n-gram set for each sequences in sequences
  sequences -- list of sequences
  n -- n-gram size"""
  for s in sequences: 
    s.calc_ngrams(n)
  return sequences # just for convenience

  
def newlines(text):
  """Unifies the new line codes"""
  text = text.replace('\r\n','\n')
  text = text.replace('\n\r','\n')
  text = text.replace('\r','\n')    
  return text
   
   
def load_fasta(filepath):
  """Loads a file with sequences in FASTA format and returns a list
  of sequences
  filepath -- path to sequence file
  """
  def parse(segment):
    """Parses a segment of a text that contains a sequence in FASTA format"""
    lines = newlines(segment).split('\n')
    name = lines[0].split()[0]
    letters = "".join(lines[1:])
    seq = Seq(name, letters.replace('*',''))
    return seq

  f = open(filepath)
  segments  = f.read().split(">")
  f.close()
  return [parse(segment) for segment in segments[1:]]


def load_phylip(filepath):  
  """Loads a file with sequences in Phylip format
  filepath -- path to Phylip file
  """
  f = open(filepath)
  lines = newlines(f.read())
  lines = lines.split('\n')[1:]
  seqs = []
  for line in lines:
    name = line[:10].strip()
    data = line[10:].replace('-','').upper()
    seq = Seq(name, name, data)
    seqs.append(seq)
  f.close()    
  return seqs    



if __name__ == "__main__":
  """Just a usage example"""
  sequences = load_sequences("data/MR_GR.ffa")
  for sequence in sequences:
    print sequence.name(), len(sequence)
  
