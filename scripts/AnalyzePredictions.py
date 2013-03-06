"""AnalyzePredictions

Usage:
  AnalyzePredictions.py --gtruth=<gtdir> --other=<odir> --parana=<pfile>

Options:
  -g, --gtruth=<gtdir>     The directory containing the ground truth networks.
  -o, --other=<odir>       The directory containing the predictions of the other method.
  -p, --parana=<pfile>     The PARANA output file.
"""
from docopt import docopt

import networkx as nx
import numpy as np
import glob
import itertools
import croc
import sys
import argparse
import evaluation
from util import allPairs
from collections import namedtuple

def graphWithCutoff( fname, cutoff=0.0 ):
    '''
    Reads in the PARANA format file 'fname' and returns a graph 
    containing all edges with weight >= cutoff.
    '''
    edges = []
    with open(fname) as ifile:
        for l in ifile:
            toks = l.rstrip().split('\t')
            edges.append( (toks[0], toks[1], float(toks[3])) )
    G = nx.Graph()
    G.add_weighted_edges_from( [ e for e in edges if e[2] > cutoff ] )
    return G

def filterGraphWithCutoff(G, cutoff=0.5):
    '''
    Returns a new graph on the same node set as 'G', but removing any edges
    with a weight < 'cutoff'.
    '''
    G2 = nx.Graph()
    G2.add_nodes_from(G.nodes())
    for u,v,d in G.edges_iter(data=True):
        if d['weight'] >= cutoff:
            G2.add_edge(u,v)
    return G2

def filteredGraph(G, nodes, cutoff=0.5):
    '''
    Returns a new graph G', on the node set 'nodes' having 
    only those edges in nodes X nodes that exist in G 
    with a weight >= cutoff.
    '''
    NG = nx.Graph()
    NG.add_nodes_from(nodes)
    for u,v in G.edges_iter():
        if ((u in nodes) and (v in nodes)) and G[u][v]['weight'] >= cutoff :
            NG.add_edge(u,v)
    return NG


def writeStatsToFile( gfname, sfname, tgraph ):
    """
    gfname is the name of the file with the predicted edgelist
    sfname is the name of the output file
    tgraph is the name of the ground truth graph
    """
    ParProbG = graphWithCutoff(gfname, 0.0)
    with open(sfname,'wb') as ofile:
        for u,v in itertools.combinations( tgraph.nodes(), 2 ):
            ofile.write("{0} {1}\n".format( ParProbG[u][v]['weight'] if ParProbG.has_edge(u,v) else 0.0, 1 if tgraph.has_edge(u,v) else 0) )

ROC_Results = namedtuple('ROC_Results', ['BEDROC', 'AUROC', 'AUPR'])
CutoffPoint = namedtuple('CutoffPoint', ['cutoff', 'tp', 'tn', 'fp', 'fn', 'precision', 'recall', 'f1'], verbose=False)

def F1Score(cutoff, tp, tn, fp, fn):
    '''
    Computes the F1-Score at the given cutoff.
    '''
    try:
        precision = tp / float(tp+fp)
        recall = tp / float(tp+fn)
        f1 = 2.0 * ((precision*recall)/(precision+recall))
        return CutoffPoint(cutoff,tp,tn,fp,fn,precision,recall,f1)
    except ZeroDivisionError as detail:
        return None

def getPerformanceAtCutoff( pgraph, tgraph, cutoff ):
    '''
    Returns a (CutoffPoint, ROC_Results) pair where the information in the cutoff
    point (e.g. the F1-Score) has been evaluated at the supplied 'cutoff'.
    '''
    import scikits.learn 
    import scipy.special
    import warnings
    import math

    sd = croc.ScoredData()
    nedges = math.ceil(scipy.special.binom(tgraph.order(),2)) + tgraph.order()
    ys = np.zeros(nedges)
    ps = np.zeros(nedges)
    for i,(u,v) in enumerate(allPairs( tgraph.nodes() )):
        y = 1 if tgraph.has_edge(u,v) else 0
        p = pgraph[u][v]['weight'] if pgraph.has_edge(u,v) else 0.0
        ys[i], ps[i] = y, p

    sd = croc.ScoredData( zip(ps,ys) )

    warnings.simplefilter("ignore")
    Ps, Rs, _ = scikits.learn.metrics.precision_recall_curve(ys, ps)
    warnings.simplefilter("always")

    sweep = list(sd.sweep_threshold(tie_mode="smooth"))
    
    vals = []
    for k,v in sd.score_labels.iteritems():
        vals = vals + [ (k,x) for x in v ]

    vals = sorted(vals, key = lambda x: x[0], reverse=True)

    l = len(vals)    
    bestCutoff = CutoffPoint(cutoff=0.0, tp=0, tn=0, fp=0, fn=0, precision=0.0, recall=0.0, f1=0.0)
    for i,(tp,tn,fp,fn) in enumerate(sweep):
        if vals[i+1][0] < cutoff:
            c = vals[i][0]
            bestCutoff = F1Score(c, tp, tn, fp, fn)
            break

    F = lambda sweep: croc.ROC(sweep).transform(croc.Linear())
    C = F(sweep)
    AUPR = scikits.learn.metrics.auc(Ps, Rs)
    return (bestCutoff, ROC_Results( BEDROC=croc.BEDROC(sd, 20.0)['BEDROC'], 
                                     AUROC=C.area(), AUPR=AUPR ))


def getOptimalCutoff( pgraph, tgraph ):
    '''
    Return a CutoffPoint that contains the optimal cutoff for this predictor
    (i.e. the cutoff at which the F1-Score is the highest).
    '''
    sd = croc.ScoredData()
    for u,v in allPairs( tgraph.nodes() ):
        sd.add( pgraph[u][v]['weight'] if pgraph.has_edge(u,v) else 0.0, 1 if tgraph.has_edge(u,v) else 0 )

    sweep = sd.sweep_threshold(tie_mode="smooth")
    vals = []
    for k,v in sd.score_labels.iteritems():
        vals = vals + [ (k,x) for x in v ]
    vals = sorted(vals, key = lambda x: x[0], reverse=True)

    def F1Score(cutoff, tp, tn, fp, fn):
        try:
            precision = tp / float(tp+fp)
            recall = tp / float(tp+fn)
            f1 = 2.0 * ((precision*recall)/(precision+recall))
            return CutoffPoint(cutoff,tp,tn,fp,fn,precision,recall,f1)
        except ZeroDivisionError as detail:
            return None

    l = len(vals)
    bestCutoff = CutoffPoint(cutoff=0.0, tp=0, tn=0, fp=0, fn=0, precision=0.0, recall=0.0, f1=0.0)
    for i,(tp,tn,fp,fn) in enumerate(sweep):
        if i < l:
            c = vals[i][0]            
            currentCutoff = F1Score(c, tp, tn, fp, fn)
            if (currentCutoff is None): continue
            if (currentCutoff.f1 >= bestCutoff.f1 or bestCutoff.cutoff == 0.0):
                bestCutoff = currentCutoff

    return bestCutoff.cutoff, bestCutoff.f1

def getCurve( pgraph, tgraph ):
    '''
    Gets the BEDROC and AUROC scores for the given prediction.
    '''
    F = lambda sweep: croc.ROC(sweep).transform(croc.Linear())
    sd = croc.ScoredData()
    for u,v in allPairs( tgraph.nodes() ):
        sd.add( pgraph[u][v]['weight'] if pgraph.has_edge(u,v) else 0.0, 1 if tgraph.has_edge(u,v) else 0 )

    sweep = sd.sweep_threshold(tie_mode="smooth")
    C = F(sweep)

    return ROC_Results( BEDROC=croc.BEDROC(sd, 20.0)['BEDROC'], AUROC=C.area() )

def findSuggestedCutoff( I, cut ):
    '''
    Returns the cutoff of edge weights on the graph 'I' at which
    'cut' fraction of all the edge weight mass is included.  For example,
    if cut is 0.9 and S is the sum of edge weights of I, this functio will
    return the smallest cutoff of edge weights where the graph generated
    by this cutoff contains a mass of at least 0.9 * S.
    '''
    #SI = I.subgraph(G.nodes())
    sw = sorted( [d['weight'] for u,v,d in I.edges_iter(data=True)], reverse=True )
    s = sum(sw)
    sums = np.array([ sum( sw[:i] ) / s for i in xrange(1,len(sw)) ])
    cut = min(cut, sums[-2])
    return sw[ np.where(sums > cut)[0][0] ]

def getMassAtCutoff( G, c ):
    sw = sorted( [d['weight'] for u,v,d in G.edges_iter(data=True)], reverse=True )
    s = sum(sw)
    cutoffSum = sum( [w for w in sw if w > c] )
    return cutoffSum / s

def cutoffForDensity(I, density):
    sw = sorted( [d['weight'] for u,v,d in I.edges_iter(data=True)], reverse=True )
    nnodes = float(I.order())
    numedges = int(density * nnodes)
    print(numedges)
    return sw[ numedges ]


def printPerf(methodName, perf, roc):
    '''
    Print the performance statistics for the given method.
    '''
    print("{:-^40}".format("Method : "+methodName))
    print("Precision = {0}, Recall = {1}, F1-Score = {2}".format(perf.precision, perf.recall, perf.f1))
    print("BEDROC: {0}, AUROC: {1}, AUPR: {2}".format(roc.BEDROC, roc.AUROC, roc.AUPR))


def main(options):
    print(options)
    gtruth = options['--gtruth']
    other = options['--other']
    parana = options['--parana']

    # The ground truth graphs (obtained from sequence similarity)
    gtGraphNames = glob.glob("{0}/*.prob".format(gtruth))
    
    """
    Dictionary mapping each species name to it's ground truth graph
    The ground truth graph is obtained by filtering the ".prob" graph 
    which contains a non-zero probability for each potential edge
    with a cutoff of 0.5.  This ensures we keep all of the nodes, even if they have no
    incident edges
    """
    gtGraphs = { fn.split("/")[-1][:-5] : \
                 filterGraphWithCutoff(nx.read_weighted_edgelist(fn)) for fn in gtGraphNames }

    """
    The graphs that are the output of the other program (in this case Pinney et al.).
    They have the same name as the ground truth files above except '.prob' is replaced
    with '.out.ppi'.
    """
    oGraphNames = [ "{0}/{1}.out.ppi".format(other, k) for k in gtGraphs.keys() ]
    oGraphs = { fn.split("/")[-1][:-8] : nx.read_weighted_edgelist(fn) for fn in oGraphNames }

    inputGraph = nx.read_adjlist('../../Parana2Data/bZIPData/input/bZIP_n0.adj')

    # Read in the Parana graph from the file    
    paranaGraph = graphWithCutoff(parana, 0.0)

    #density = inputGraph.size() / float(inputGraph.order())

    # inputGraphNames = glob.glob("{0}/bZIP*input.adj".format(options.other))
    # inputGraph = nx.read_adjlist(inputGraphNames[0])
    # print("input graph order = {0}, size = {1}".format(inputGraph.order(), inputGraph.size()))
    # cutoff = 0.99
    # c = findSuggestedCutoff( paranaGraph, inputGraph, cutoff )
    # cOpt = getOptimalCutoff( paranaGraph, inputGraph )
    # print("Optimal cutoff = {0}".format(cOpt))
    # evaluation.printStats( filteredGraph(paranaGraph, inputGraph.nodes(), cutoff=c ), inputGraph )
    # print >>sys.stderr, "Parana 2.0    : {0}".format(getCurve(paranaGraph, inputGraph))

    for i, (gtName, gtGraph) in enumerate(gtGraphs.iteritems()):
        cutoff = 0.99

        pGraph = paranaGraph.subgraph(gtGraph.nodes())
        #c = findSuggestedCutoff(pGraph, cutoff ) 
        #c = cutoffForDensity(pGraph, density)
        
        # Find the best cutoff for the parana graph rather than 
        # using the heuristic above
        c1, bf1 = getOptimalCutoff( oGraphs[gtName], gtGraph )      
        c, bf1 = getOptimalCutoff( paranaGraph, gtGraph )      
        #print("optimal cutoff is {0}".format(c))
        #c = 0.05
        #c  = 0.03

        specString = "Species : {0}".format(gtName)
        print("{:=^80}".format(specString))
        
        (otherPerf, otherROC) = getPerformanceAtCutoff( oGraphs[gtName], gtGraph, c1 ) #0.5 )
        (paranaPerf, paranaROC) = getPerformanceAtCutoff( pGraph, gtGraph, c )
        
        printPerf("Pinney et al.", otherPerf, otherROC)
        printPerf("Parana 2.0", paranaPerf, paranaROC)

        # Old evaluation code
        #evaluation.printStats( filteredGraph(oGraphs[gtName], gtGraph.nodes()), gtGraph )
        #print >>sys.stderr, "Pinney et. al : {0}".format(getCurve(oGraphs[gtName], gtGraph))
        
        # Old evaluation code
        #evaluation.printStats( filteredGraph(paranaGraph, gtGraph.nodes(), cutoff=c ), gtGraph )
        #print >>sys.stderr, "Parana 2.0    : {0}".format(getCurve(paranaGraph, gtGraph))

        print("\n")
    sys.exit(0)

if __name__ == "__main__" :
    arguments = docopt(__doc__, version="AnalyzePredictions 1.0")
    main(arguments)