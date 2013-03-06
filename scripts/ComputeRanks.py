"""
Usage: ComputeRanks.py PDIR OGROUP (heatmap|cdf|pdf) [--relative] [--noself] -o FILE

Arguments:
  PDIR         Directory containing cross-validation predictions
  OGROUP       File containing orthology group mapping
Options:
  -h --help    show this
  --relative   compute relative (rather than absolute) ranks
  --noself     don't test self-loops
  -o FILE      specify the output file [default: plot.pdf]
"""
from docopt import docopt

import os
import sys
import glob
import itertools
from GenMultipleAlignments import parseOrthoFile
import cogent
import subprocess
import numpy as np
import networkx as nx
from subprocess import call
from collections import namedtuple
import pylab
from docopt import docopt
from ete2 import Phyloxml

ScoredEdge = namedtuple('ScoredEdge', ['p1', 'p2', 'score'])

extantName = "../../Parana2Data/HerpesPPIs/extant.adj"
#extantName = HerpesExtant.adj"

temp = nx.read_adjlist(extantName)
nmap = { n : n.upper() for n in temp.nodes_iter() }
ExtantNetwork = nx.Graph()
ExtantNetwork.add_nodes_from( [nmap[n] for n in temp.nodes_iter()] )
ExtantNetwork.add_edges_from( [(nmap[u], nmap[v]) for u,v in temp.edges_iter()] )

def combinationsWithSelf(l):
  return list(itertools.combinations(l,2)) + [ (x,x) for x in l ]

def readScoreFile(fname, noself, randomize=False):

  # The strings naming the proteins whose interaction was removed in
  # this input
  tstring = fname.split('@')[-2].split("#")
  # Convert to upper case and make into an edge name
  tedge = (tstring[0].upper(), tstring[1].upper())

  # Read in the phylogenies for the orthology groups
  treeDir = "../../Parana2Data/HerpesPPIs/trees/rearranged"#"dataOut_June17/rearranged"
  baseFile = fname.split('/')[-1]
  orthoGroup1 = baseFile.split('@')[0]
  orthoGroup2 = baseFile.split('@')[1]
  t1 = "{0}/{1}.xml.rooting.0.ntg.reconciled.0.ntg.rearrange.0.ntg".format(treeDir,orthoGroup1)
  t2 = "{0}/{1}.xml.rooting.0.ntg.reconciled.0.ntg.rearrange.0.ntg".format(treeDir,orthoGroup2)

  if ( not (os.path.exists(t1) and os.path.exists(t2)) ):
    return None, None

  # The extant (non ancestral, non lost) nodes from the two homology groups
  getName = lambda x : x.Name.upper()
  pxml = Phyloxml()
  pxml.build_from_file(t1)
  pxml.build_from_file(t2)
  la = filter( lambda x : x.find('LOST') == -1, map( lambda x: x.upper(), pxml.phylogeny[0].get_leaf_names() ) )
  lb = filter( lambda x : x.find('LOST') == -1, map( lambda x: x.upper(), pxml.phylogeny[1].get_leaf_names() ) )

  # The set of all possible interactions among the two homology groups
  #pe = list(itertools.product(la,la)) + list(itertools.product(la,lb)) + list(itertools.product(lb,lb))
  possibleEndpoints =  combinationsWithSelf(la) + list(itertools.product(la,lb)) + combinationsWithSelf(lb) \
                       if not (orthoGroup1 == orthoGroup2 ) else combinationsWithSelf(la)

  # From among all possible endpoints, only those protein pairs that reside in the
  # same species represent a potential edge                      
  allPossibleEdges = filter( lambda (x,y): x.split('_')[-1] == y.split('_')[-1], possibleEndpoints )
    
  # an edge has both endpoints in the set of extant nodes
  inCurrentGroups = lambda e, x, y: (e[0] in x or e[0] in y) and (e[1] in x or e[1] in y)
  # an edge is relevant if it's constrained to the current groups
  relevantExtantEdges = [ e for e in ExtantNetwork.edges_iter() if inCurrentGroups(e,la,lb) ]
  # the set of potential edges that don't appear in the input network
  nonPresentEdgesMinusTarget = list(set([ x for x in allPossibleEdges if not ExtantNetwork.has_edge(x[0],x[1])]))
  # the same as above but including our target edge
  nonPresentEdges = nonPresentEdgesMinusTarget + [tedge]

  import random

  # Ancestral edges start with an N or R
  ancestral = ['R','N']

  # The node is valid if it is neither lost nor ancestral
  isValidNode = lambda x : (x[0] not in ancestral) and (x.find('LOST') == -1)

  # Is the edge u,v the target edge?
  isCurrentEdge = lambda u,v : (u == tedge[0] and v == tedge[1]) or (u == tedge[1] and v == tedge[0])
  isRealEdge = lambda u,v : (not isCurrentEdge(u,v)) and ExtantNetwork.has_edge(u,v)
  isValidEdge = lambda u,v : ((isValidNode(u) and isValidNode(v)) and (not isRealEdge(u,v)))

  # Is u,v one of the edges we wish to consider?
  def inPotentialEdges(u,v) :
    contains = (u,v) in nonPresentEdges or (v,u) in nonPresentEdges
    if noself:
      return u != v and contains
    else:
      return contains

  def isEdge( se, p1, p2 ):
      r =  ((se.p1 == p1 and se.p2 == p2) or (se.p1 == p2 and se.p2 == p1))
      return r

  scoredEdges = []

  nonEdgesWithProb = set( nonPresentEdges )
  with open(fname,'rb') as ifile:
    for l in ifile:
      toks = l.rstrip().split()
      p1 = toks[0].upper()
      p2 = toks[1].upper()
      s = float(toks[3])
      if inPotentialEdges(p1,p2):
        #s = random.uniform(0.0,1.0)
        se = ScoredEdge(p1,p2,s)
        scoredEdges.append( se )
        nonEdgesWithProb.discard((p1,p2))
        nonEdgesWithProb.discard((p2,p1))

  rev = True
  for u,v in (nonEdgesWithProb - set(nonPresentEdges)):
    scoredEdges.append(ScoredEdge(u, v, 0.0))
  
  # cost = 0.0
  # for u,v in nonPresentEdges:
  #     se = ScoredEdge(u,v,cost)
  #     fe = [ e for e in scoredEdges if isEdge(e, u, v) ]
  #     if len(fe) == 0:
  #         scoredEdges.append(se)

  random.shuffle(scoredEdges)
  scoredEdges = list(enumerate(sorted( scoredEdges, key=lambda x: x.score, reverse=rev )))
  print(len(scoredEdges))
  print(t1,t2)
  print("Target Edge = {0}".format(tedge))
  print("Extant Edges = {0}".format(relevantExtantEdges))
  print("Potential Edges = {0}".format(nonPresentEdges))
  print("Scored Edges = {0}".format(scoredEdges))

  res = [ x for x in scoredEdges if isEdge(x[1], tedge[0], tedge[1])  ]

  if len(res) > 0:
    print(res)
    # Prev (ISMB)
    #print(res[0][0],float(len(nonPresentEdges)-1))
    #return (res[0][0], float(len(nonPresentEdges)-1))
    # New
    print(res[0][0],float(len(scoredEdges)-1))
    return (res[0][0], float(len(scoredEdges)-1))
    
  else:
    raise 'Hell'
    #rr = random.uniform(len(res)/float(len(nonPresentEdges)+1),1.0)
    #return (rr,1.0)

def pairwiseHeatMap(rankDict, orthoInds):
  import numpy as np
  import pylab
  from matplotlib import cm as CM

  hmMap = np.zeros( (len(orthoInds), len(orthoInds)) )
  hmMap[:,:] = np.nan

  for k,v in rankDict.iteritems():
    gi, gj = k
    hmMap[orthoInds[gi]][orthoInds[gj]] = hmMap[orthoInds[gj]][orthoInds[gi]] = np.mean(v)

  pylab.imshow(hmMap, cmap=CM.jet,  interpolation='nearest')
  #heatmap[:, ::-1], origin='lower', extent=extent, aspect=2, interpolation='nearest')
  pylab.clim(0,1.0)
  pylab.colorbar()
  print(rankDict)
  pylab.ylabel('Orthology Groups', labelpad=25)
  pylab.xlabel('Orthology Groups', labelpad=25)
  #pylab.hexbin(locs, counts, gridsize=20, cmap=CM.jet)

def heatMap( rankDict ):
  import numpy as np
  import pylab
  from matplotlib import cm as CM

  locs = []
  counts = []
  for i,(k,v) in enumerate(rankDict.iteritems()):
    locs.append(i)
    counts.append( np.mean(v) )

  heatmap, xedges, yedges = np.histogram2d(locs, counts, bins=20)
  heatmap = np.rot90(heatmap,3)
  extent = [xedges[1], xedges[-1], yedges[1], yedges[-1]]# max(yedges)]
  print(extent)
  pylab.xticks([])
  pylab.imshow(heatmap[:, ::-1], origin='lower', extent=extent, aspect=2, interpolation='nearest')
  pylab.ylabel('Absolute Ranks')
  pylab.xlabel('Pairs of Orthology Groups')
  pylab.title('Average Ranks Across Orthology Group Pairs')
  #pylab.hexbin(locs, counts, gridsize=20, cmap=CM.jet)

def boxPlot( rankDict ):
  import numpy as np

  pylab.rcParams.update({'lines.linewidth' : 2.0})

  counts = []
  for i,(k,v) in enumerate(rankDict.iteritems()):
    counts.append(v)
  pylab.boxplot(counts)
  pylab.xticks(np.arange(1,len(rankDict)+1),rankDict.keys(), rotation=60)
  #pylab.ylim(0.0, 1.0)
  pylab.title('Rank Distributions by Orthology Group')
  pylab.xlabel('Orthology Classes')
  pylab.ylabel('Relative Ranks')

def main():
  import numpy
  import scipy

  arguments = docopt(__doc__, version='ComputeRanks v1.0')
  relative = arguments['--relative']
  print(arguments)

  params = { 'font.size' : 20,
             'lines.linewidth' : 3.0,
  }

  pylab.rcParams.update(params)

  import GenMultipleAlignments
  orthoDict = GenMultipleAlignments.parseOrthoFile(arguments['OGROUP'])[1]
  revOrthoDict = {}

  numOrthoGroups = len(orthoDict)
  orthoInds = { k:i for i,k in enumerate(orthoDict.iterkeys()) }

  for orthoGroup, spDict in orthoDict.iteritems():
    for spName, geneList in spDict.iteritems():
      newGeneList = [e.upper() for e in geneList]
      orthoDict[orthoGroup][spName] = newGeneList
      for e in newGeneList:
        revOrthoDict[e] = orthoGroup

  print(revOrthoDict)

  G = nx.read_adjlist(extantName)
  N = G.order()
  print(N)
  numPossible = float(N*N)
  noself = arguments["--noself"]

  tested = set()
  rankDict = {}
  ranks = []
  for fname in glob.glob("{0}/*@txt".format(arguments['PDIR'])):
    tstring = fname.split('@')[-2].split("#")
    tedge = (tstring[0].upper(), tstring[1].upper())
    
    print("TEDGE ",tedge)
    if noself and tedge[0] == tedge[1]:#tedge in tested: 
#      exit(1)
      continue
    else:
      tested.add(tedge)
      tested.add((tedge[1],tedge[0]))

    pref = tstring[0].split('_')[-1]
    
    rank, listSize = readScoreFile(fname, noself, randomize=True)
    if rank is None:
      continue

    rank = rank / float(listSize) if relative else rank
    ranks.append(rank)

    if not ((tedge[0] in revOrthoDict) and (tedge[1] in revOrthoDict)):
      continue
    print( revOrthoDict )
    print([revOrthoDict[tedge[0]], revOrthoDict[tedge[1]] ])
    pref = tuple(sorted([revOrthoDict[tedge[0]], revOrthoDict[tedge[1]]]))

    #hmMap[orthoInds[pref[0]]][orthoInds[pref[1]]] = rank
    #hmMap[orthoInds[pref[1]]][orthoInds[pref[0]]] = rank

    #for pref in [revOrthoDict[tedge[0]], revOrthoDict[tedge[1]]]:
    if pref in rankDict:
      rankDict[pref].append(rank)
    else:
      rankDict[pref] = [rank]
  
  from collections import Counter
  cnt = Counter(ranks)
  print(sorted( [(k,v) for k,v in cnt.iteritems()], key = lambda x: x[0] ))
  print( sum( [v for k,v in cnt.iteritems()] ) )
  print(rankDict)
  print("NUMBER OF EXPERIMENTS = {0}".format(len(ranks)))
  print("MEAN RELATIVE RANK = {0}".format(np.array(ranks).mean()))
  print("MEDIAN RELATIVE RANK = {0}".format(np.median(ranks)))

  if arguments['heatmap']:
    pairwiseHeatMap(rankDict, orthoInds)

  elif arguments['cdf']:
    print len(rankDict)
    print rankDict
    a = numpy.array(ranks)
    print(ranks)
    num_bins = len(a)
    counts, bin_edges = numpy.histogram(a, bins=num_bins, normed=False)
    print("bin edges",bin_edges)
    cdf = numpy.cumsum(counts, axis=0)
    bin_edges = list(bin_edges) #+ [1.0]
    cdf = list(cdf) #+ [cdf[-1]]
    # relative
    if relative:
        tot = numpy.max(cdf)
        cdf = [ c / float(tot) for c in cdf ]
        pylab.plot(bin_edges[1:] + [1.0], cdf + [cdf[-1]] )
        pylab.xlabel('Relative Rank of Imputed Edge')
        pylab.ylabel('Fraction of Edges With Relative Rank < x')
        pylab.grid(axis='both')
        pylab.xlim((0.0,1.0))
        pylab.ylim((0.0,1.0))
        print(cdf)
        print(bin_edges)
    else:
      # absolute
      pylab.plot(bin_edges[1:], cdf )

  elif arguments['pdf']:
      a = numpy.array(ranks)
      num_bins = 5
      # relative
      if relative:
          print(a)
          pylab.hist(a,bins=num_bins)
          #pylab.xticks( [0.1] + [1.3, 1.0/2.0, 0.7] + [0.9] , ('0', '1/3', '1/2', '2/3', '1') )
          pylab.xlim((0.0,1.0))
          pylab.xlabel('Relative Rank of Imputed Edge')
          pylab.ylabel('# of Imputed Edges with Given Relative Rank', labelpad=25)
          pylab.grid(axis='both')
          #pylab.xlim((0.0,1.0))
          #pylab.ylim((0.0,1.0))
  import matplotlib
  ax = pylab.gca()
  #ax.title.set_y(1.05)
  #pylab.tight_layout()
  pylab.savefig(arguments['-o'])

  #heatMap(rankDict)



  #boxPlot(rankDict)

#  yticks = pylab.yticks()[0]
#  for yt in yticks:
#    pylab.axhline(yt, linewidth=1.0, color='gray', linestyle='dashed')

#  xticks = pylab.xticks()[0]
#  for xt in xticks:
#    pylab.axvline(xt, linewidth=1.0, color='gray', linestyle='dashed')

#  pylab.hist(a, bins=num_bins, normed=False)
#  pylab.title("Distribution of ranked list lengths")
#  pylab.xlabel("Length of Edge List")
#  pylab.ylabel("# of Lists With Given Length")
#  pylab.xlim(bin_edges[0], bin_edges[-1])

#  pylab.plot(relRanks)
#  pylab.ylim(0, len(a))
#  pylab.xlim(0.0,1.0)


#  pylab.title("Cumulative Count of Ranks)
#  pylab.ylabel("# of Experiments")
#  pylab.xlabel("Absolute Rank")

#  pylab.scatter( [x[0] for x in relRanks], [x[1] for x in relRanks] )
#  pylab.xlim(0, max([x[1] for x in relRanks]))
#  pylab.ylim(0, max([x[1] for x in relRanks]))
#  pylab.title("Rank vs. List Length")
#  pylab.xlabel("Rank")
#  pylab.ylabel("Rank List Length")
  pylab.show()

if __name__ == "__main__":
  main()
