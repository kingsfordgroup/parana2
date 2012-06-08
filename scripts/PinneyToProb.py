import networkx as nx
import numpy as np
import itertools
import argparse
import datetime
import scipy.stats

def filteredGraph( G, cutoff=0.0 ):
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for u,v,d in G.edges_iter(data=True):
        if d['weight'] > cutoff:
            H.add_edge(u,v,weight=d['weight'])
    return H

def createParser():
    parser = argparse.ArgumentParser(prog="./AnalyzePredictions")
    parser.add_argument("-i", "--input", default=None, help="Pinney input data")
    parser.add_argument("-o", "--output", default=None, help="Output probabilistic data")
    parser.add_argument("-l", "--logistic", default=None, help="Convert scores to probabilities using a logistic function")
    parser.add_argument("-b", "--binary", action='store_true' , help="binarize the data")
    return parser

def writeAdjList(G, adjDict, lp, ofname):

    with open(ofname,'wb') as ofile:
        dt = datetime.datetime.now()
        ofile.write("# File written on: {0}\n".format(dt.strftime("%A, %d. %B %Y %I:%M%p")))
        ofile.write("# Logistic parameter: {0}\n".format(lp))
        # write out each node
        for n in G.nodes_iter():
            # if it has neighbors
            if n in adjDict:
                nlist = list(adjDict[n])
                numNeighbors = len(nlist)
                ofile.write("{0} {1}\n".format(n, numNeighbors))
                for l in nlist:
                    ofile.write("{0} {1}\n".format(l[1], l[2]))
            else:
                ofile.write("{0} 0\n".format(n))

    return True

def main():
    parser = createParser()
    options = parser.parse_args()


    G = nx.read_weighted_edgelist(options.input, delimiter='\t', comments='#')
    weights = [ d['weight'] for u,v,d in G.edges_iter(data=True) ]
    minW = min(weights)
    maxW = max(weights)
    cutoffW = 30.6 # (30.6 - minW) / (maxW - minW)
    print(minW, maxW, cutoffW)

    if options.binary :
        I = filteredGraph(G, cutoffW)
        nx.write_adjlist(I, options.output)
    else:
        import numpy

        lp = float(options.logistic)

        I = filteredGraph(G, cutoffW)

        stddev = np.std( [cutoffW - d['weight'] for u,v,d in I.edges_iter(data=True)]  + [d['weight'] - cutoffW for u,v,d in I.edges_iter(data=True)] )

        n = scipy.stats.norm( loc = cutoffW, scale = stddev)
        print(stddev)

        #def prob(x):
        #    return min(1.0, 2.0*n.cdf(x))#(0.5 + ((x-cutoffW) / (maxW-cutoffW))) if (x >= cutoffW) else (0.5*(x-minW) / (cutoffW-minW))

        prob = lambda x : 1.0 / ( 1.0 + np.exp( lp * (-x + cutoffW) ) )

        H = nx.Graph()
        H.add_nodes_from(G.nodes())

        for u,v,d in G.edges_iter(data=True):
            w = prob(d['weight'])
            H.add_edge( u, v, weight=w )

        strs = sorted([ (u,v,d['weight']) for u,v,d in H.edges_iter(data=True) ])
        adjDict = { k : list(v) for k,v in itertools.groupby(strs, lambda x : x[0] ) }
        writeAdjList(H, adjDict, lp, options.output)

if __name__ == "__main__":
    main()