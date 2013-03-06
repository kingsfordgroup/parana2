"""Usage: BuildSingleTree.py <network> <seqfile> -o <ofile>

-h --help    show this
-o ofile     specify the output file [default: all.nwk]
"""

import sys
import os
import glob
import errno

import cogent
import networkx as nx
from cogent.app.muscle import align_unaligned_seqs
from cogent.app.fasttree import build_tree_from_alignment
from docopt import docopt

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def main():
    arguments = docopt(__doc__, version='BuildSingleTree v1.0')
    print(arguments)
    netname = arguments['<network>']
    seqname = arguments['<seqfile>']
    ofname = arguments['-o']

    print("creating tree for {0}".format(netname))
    print("using sequences from {0}".format(seqname))

    G = nx.read_adjlist(netname)
    try:
        seqs = cogent.LoadSeqs(seqname, moltype=cogent.PROTEIN, aligned=False)
    except Exception as e:
        print(e)
        sys.exit(0)

    aln = align_unaligned_seqs(seqs, cogent.PROTEIN)
    t = build_tree_From_alignment(aln, cogent.PROTEIN)
    print("tree  = {0}".format( str(t) ) )
    with open( os.path.sep.join( [ "all" ] ), 'wb') as ofile:
        ofile.write( str(t).replace("'","") )

      # fdir = sys.argv[1]
      # odir = sys.argv[2]
      # mkdir_p(odir)

      # for fname in glob.iglob("{0}/*.fasta".format(fdir)):

      #   groupName = fname.split( os.path.sep )[-1].rstrip(".fasta")
      #   print("group {0}".format(groupName))
      #   try:
      #     seqs = cogent.LoadSeqs(fname, moltype=cogent.PROTEIN, aligned=False)
      #   except Exception as e:
      #     print(e)
      #     exit(0)
      #   aln = align_unaligned_seqs(seqs, cogent.PROTEIN)
      #   t = build_tree_from_alignment(aln, cogent.PROTEIN)
    #   print("tree for group {0}".format( str(t) ) )
    #   with open( os.path.sep.join( [ odir, groupName+".nwk" ] ), 'wb') as ofile:
    #       ofile.write( str(t).replace("'","") )

if __name__ == "__main__": main()
