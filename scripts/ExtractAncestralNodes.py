import os
import sys
import glob
from GenMultipleAlignments import parseOrthoFile
import cogent

def main():
  (spec, odict) = parseOrthoFile( sys.argv[1] ) 
  recDir = sys.argv[2]
  ofname = sys.argv[3]
  
  orthoGroups = odict.keys()
  with open(ofname,'wb') as ofile:
    for og in orthoGroups:
      ifname = "{0}/{1}.nwk.rooting.0".format(recDir,og)
      t = cogent.LoadTree(ifname)
      ofile.write(t.Name+"\n") 

  sys.exit(1)

if __name__ == "__main__" : 
  main()
