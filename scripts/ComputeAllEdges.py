"""
Usage: ComputeAllEdges.py --method=<m> --recdir=<r> --ortho=<o> --outdir=<od> [--beta=<b>]

Options:
--method=<m>  Method to use to compute the cross validation results {Parana | RWS | CV}
--recdir=<r>  Folder containing reconciled trees
--outdir=<o>  Folder to which results should be written
--ortho=<o>   Orthogy group file
--beta=<b>    Beta for Parana method [default: 80.0]
"""
from docopt import docopt
import os
import sys
import glob
import itertools
from GenMultipleAlignments import parseOrthoFile
import cogent
import subprocess
import networkx as nx
import numpy as np
import tempfile
from subprocess import call
import os, errno
from ete2 import Phyloxml
from progressbar import ProgressBar
from lxml import etree as ET

sys.path.append("../../rws")

extantName = "../../Parana2Data/HerpesPPIs/extant_restricted.adj"

def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc: # Python >2.5
    if exc.errno == errno.EEXIST:
      pass
    else: raise

def iterWithSelf(groups):
    return itertools.chain(itertools.combinations(groups,2), [(g,g) for g in groups])

def constructParanaCall(gfile, ofilename, beta, og1, og2, t1name, t2name):
  parana = "../build/bin/parana2"
  callargs = [parana, "-t {0}".format(gfile.name), \
  "pars", \
  "-u", "-o", ofilename,\
  "-k", "40", "-r", "1.2","-p", "1.0", "-b", str(beta)]

  # Add the first orthology group
  callargs += ["-d", t1name]

  # If the two orthology groups are different, add the second
  if og1 != og2 :
    callargs += ["-d", t2name]

  return callargs

def runRWS(G, ofilename):
  import rwsPredict
  vertToInd = { v:i for i,v in enumerate(G.nodes()) }
  A = np.array(nx.adjacency_matrix(G))

  #A -= np.diag(np.diag(A))
  #A += np.eye(A.shape[0])

  # u = vertToInd[redge[0]]
  # v = vertToInd[redge[1]]
  # A[u][v] = 0.0
  # A[v][u] = 0.0

  G2 = rwsPredict.corrGraph(G, A, filterSelf=False)
  with open(ofilename,'wb') as ofile:
    rwsPredict.writeInParanaFormat(G2, ofile)

def main(args):
  usingParana = (args["--method"].upper() == "PARANA")
  usingRWS = (args["--method"].upper() == "RWS")
  outputCVFile = (args["--method"].upper() == "CV")
  (spec, odict) = parseOrthoFile( args["--ortho"] ) #sys.argv[1] )
  recDir = args["--recdir"]#sys.argv[2]
  ofname = args["--outdir"]#sys.argv[3]
  beta = float(args["--beta"]) 
  
  if not outputCVFile:
    mkdir_p(ofname)

  getName = lambda x : x.Name

  def getRelevantEdges( adjGraph, t1, t2 ):
    pxml = Phyloxml()
    pxml.build_from_file(t1)
    pxml.build_from_file(t2)
    la = filter( lambda x : x.find('LOST') == -1, pxml.phylogeny[0].get_leaf_names() )
    lb = filter( lambda x : x.find('LOST') == -1, pxml.phylogeny[1].get_leaf_names() )
    #lb = filter( lambda x : x.find('LOST') == -1, map( getName, cogent.LoadTree(t2).tips() ) )
    crossValidationEdges = filter( lambda (x,y) : ((x in la) and (y in lb)) or ((y  in la) and (x in lb))  , adjGraph.edges() )
    relevantEdges = filter( lambda (x,y) : ((x in la) or (x in lb)) and ((y in la) or (y in lb)) , adjGraph.edges() )
    newGraph = nx.Graph()
    newGraph.add_nodes_from( la + lb )
    newGraph.add_edges_from( relevantEdges )
    return newGraph, crossValidationEdges


  if outputCVFile:
    root = ET.Element("cvtest", name="herpes_loocv")
    
  ctr = 0
  orthoGroups = odict.keys()
  progress = ProgressBar(len(orthoGroups)*len(orthoGroups)).start()
  for i, (og1, og2) in enumerate(iterWithSelf(orthoGroups)):
      progress.update(i)
      #print("computing ancestral state between groups {0} and {1}".format(og1, og2))
      #print("i = {0}".format(i))
      alist = nx.read_adjlist(extantName)
      

      t1name = "{0}/{1}.xml.rooting.0.ntg.reconciled.0.ntg.rearrange.0.ntg".format(recDir,og1)
      t2name = "{0}/{1}.xml.rooting.0.ntg.reconciled.0.ntg.rearrange.0.ntg".format(recDir,og2)

      relevantGraph, crossValidationEdges = getRelevantEdges(alist, t1name, t2name)
      if len(crossValidationEdges) > 1:
          for e in crossValidationEdges:
            cGraph = relevantGraph.copy()
            cGraph.remove_edge(e[0], e[1])
            ofileName = "{0}/{1}@{2}@removed@{3}#{4}@txt".format(ofname, og1, og2, e[0], e[1])
            assert(cGraph.size() == relevantGraph.size()-1)
              
            callargs = None
            if usingParana:
              with tempfile.NamedTemporaryFile(mode='wb', suffix='.adj', prefix='tmp', dir="/tmp", delete=True) as gfile:
                nx.write_adjlist(cGraph, gfile.name)
                callargs = constructParanaCall(gfile, ofileName, beta, og1, og2, t1name, t2name)
                os.system(" ".join(callargs)+"> /dev/null 2>&1")
            elif usingRWS:
              copy = alist.copy()
              copy.remove_edge(e[0],e[1])
              #print(cGraph.size())
              runRWS(copy, ofileName)
            elif outputCVFile:
              tset = ET.SubElement(root, "testset", name="fold_{0}".format(ctr))
              u,v = sorted((e[0],e[1]))
              ET.SubElement(tset, "edge", u=u, v=v)      
              ctr += 1
            else:
              raise "Method must be one of {Parana|RWS}"

        #p = subprocess.Popen( callargs , shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE )
        #out, err = p.communicate()
      if outputCVFile:
        with open(ofname,'wb') as ofile:
          ofile.write(ET.tostring(root, pretty_print=True))
          
  progress.finish()
  sys.exit(1)

if __name__ == "__main__" :
  arguments = docopt(__doc__, version='ComputeAllEdges v1.0')
  main(arguments)
