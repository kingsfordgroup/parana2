import networkx as nx
import glob
import itertools
import croc
import sys
import argparse
import evaluation
from collections import namedtuple

def graphWithCutoff( fname, cutoff=0.0 ):
    edges = []
    with open(fname) as ifile:
        for l in ifile:
            toks = l.rstrip().split('\t')
            edges.append( (toks[0], toks[1], float(toks[3])) )
    G = nx.Graph()
    G.add_weighted_edges_from( [ e for e in edges if e[2] > cutoff ] )
    return G

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

ROC_Results = namedtuple('ROC_Results', 'BEDROC, AUROC')

def getCurve( pgraph, tgraph ):
    F = lambda sweep: croc.ROC(sweep).transform(croc.Linear())#Exponential(7.0))
    sd = croc.ScoredData()
    for u,v in itertools.combinations( tgraph.nodes(), 2 ):
        sd.add( pgraph[u][v]['weight'] if pgraph.has_edge(u,v) else 0.0, 1 if tgraph.has_edge(u,v) else 0 )

    C = F(sd.sweep_threshold(tie_mode="smooth"))
    return ROC_Results( BEDROC=croc.BEDROC(sd, 20.0)['BEDROC'], AUROC=C.area() )

def createParser():
    parser = argparse.ArgumentParser(prog="./AnalyzePredictions")
    parser.add_argument("-g", "--gtruth", default=None, help="Directory containing ground truth networks")
    parser.add_argument("-o", "--other", default=None, help="Directory containing other predicitons")
    parser.add_argument("-p", "--parana", default=None, help="File containing PARANA predictions")
    return parser

def filteredGraph(G, nodes, cutoff=0.5):
    NG = nx.Graph()
    for u,v in G.edges_iter():
        if ((u in nodes) and (v in nodes)) and G[u][v]['weight'] > cutoff :
            NG.add_edge(u,v)
    return NG

def main():
    parser = createParser()
    options = parser.parse_args()

    gtGraphNames = glob.glob("{0}/*.sim.cut".format(options.gtruth))
    gtGraphs = { fn.split("/")[-1][:-8] : nx.read_edgelist(fn) for fn in gtGraphNames }
    print(gtGraphs)

    oGraphNames = [ "{0}/{1}.out.ppi".format(options.other, k) for k in gtGraphs.keys() ]
    oGraphs = { fn.split("/")[-1][:-8] : nx.read_weighted_edgelist(fn) for fn in oGraphNames }
    print(oGraphNames)

    paranaGraph = graphWithCutoff(options.parana, 0.0)

    for gtName, gtGraph in gtGraphs.iteritems():
        print(gtName)
        print("==================")
        evaluation.printStats( filteredGraph(oGraphs[gtName], gtGraph.nodes()), gtGraph )
        print >>sys.stderr, "Pinney et. al : {0}".format(getCurve(oGraphs[gtName], gtGraph))
        evaluation.printStats( filteredGraph(paranaGraph, gtGraph.nodes(), cutoff=0.02 ), gtGraph )
        print >>sys.stderr, "Parana 2.0    : {0}".format(getCurve(paranaGraph, gtGraph))
        print("\n")
    sys.exit(0)

if __name__ == "__main__" :

    main()