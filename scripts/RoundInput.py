import networkx as nx
import itertools
import PinneyToProb
import sys

def readWeightedMultiAdj(fname):
    G = nx.Graph()
    ifile = open(fname)
    lines = filter( lambda x : not x.startswith("#"), ifile.readlines() )
    ifile.close()

    i = 0
    lenLines = len(lines)
    while i < lenLines:
        toks = lines[i].rstrip().split()
        vname, nn = toks[0], int(toks[1])
        G.add_node(vname)
        i += 1
        j = i
        mline = i+nn
        while j < mline:
            toks = lines[j].rstrip().split()
            n,w = toks[0], float(toks[1])
            j += 1
            i += 1
            G.add_edge(vname, n, weight=w)

    return G


def main():
    print(sys.argv)
    G = readWeightedMultiAdj(sys.argv[1])

    for u,v,d in G.edges_iter(data=True):
       d['weight'] = round(d['weight'],3)

    strs = sorted([ (u,v,d['weight']) for u,v,d in G.edges_iter(data=True) ])
    adjDict = { k : list(v) for k,v in itertools.groupby(strs, lambda x : x[0] ) }
    PinneyToProb.writeAdjList(G, adjDict, 0.0, sys.argv[2])


if __name__ == "__main__":
    main()