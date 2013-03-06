import sys
import os
import glob
import errno

import cogent
from cogent.app.muscle import align_unaligned_seqs
from cogent.app.fasttree import build_tree_from_alignment


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def main():
  fdir = sys.argv[1]
  odir = sys.argv[2]
  mkdir_p(odir)

  for fname in glob.iglob("{0}/*.fasta".format(fdir)):

    groupName = fname.split( os.path.sep )[-1].rstrip(".fasta")
    print("group {0}".format(groupName))
    try:
      seqs = cogent.LoadSeqs(fname, moltype=cogent.PROTEIN, aligned=False)
    except Exception as e:
      print(e)
      exit(0)
    aln = align_unaligned_seqs(seqs, cogent.PROTEIN)
    t = build_tree_from_alignment(aln, cogent.PROTEIN)
    print("tree for group {0}".format( str(t) ) )
    with open( os.path.sep.join( [ odir, groupName+".nwk" ] ), 'wb') as ofile:
        ofile.write( str(t).replace("'","") )

if __name__ == "__main__": main()
