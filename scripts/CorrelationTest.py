"""
Usage: CorrelationTest.py ENET SSCORE PRED 

Arguments:
   ENET         Extant Network
   SSCORE       Score predicted based on sequence
   PRED         Parana Predictions
Options:
   -h --help    show this
   --relative   compute relative (rather than absolute) ranks
   -o FILE      specify the output file [default: plot.pdf]
"""
from docopt import docopt

import glob
import itertools
from collections import defaultdict

import networkx as nx
import numpy as np
import scipy
import matplotlib.pyplot as plt


class Edge(object):
    def __init__(self, u, v, p):
        self._u, self._v = sorted((u,v))
        self._p = p

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self._u == other._u) and (self._v == other._v)
        else:
            False

    def __ne__(self, other): return (not self.__eq__(other))

    def __lt__(self, other):
        return (self._u, self._v) < (other._u, other._v)

    def __le__(self, other):
        return (self._u, self._v) <= (other._u, other._v)

    def __ge__(self, other):
        return (self._u, self._v) >= (other._u, other._v)

    def __str__(self):
        return "({0},{1}):{2}".format(self._u, self._v, self._p)

def main():
    arguments = docopt(__doc__, version='CorrelationTest v1.0')
    #enet = nx.read_adjlist(arguments['ENET'])
    seqfile = arguments['SSCORE']

    seqEdges = []
    with open(seqfile) as ifile:
        for l in ifile:
            toks = l.strip().split()
            if len(toks) == 3:
                u, v = sorted(toks[:2])
                p = float(toks[2])
                seqEdges.append(Edge(u, v, p))

    predEdges = []
    predfile = arguments['PRED']

    with open(predfile) as ifile:
        for l in ifile:
            u, v, f, p = l.strip().split()
            u, v = sorted([u, v])
            p = float(p)
            predEdges.append(Edge(u, v, p))

    predEdges = sorted(predEdges)
    import bisect

    import croc
    dat = croc.ScoredData()
    vals = []
    for e in seqEdges:
        i = bisect.bisect(predEdges, e)
        pedge = predEdges[i - 1]
        if pedge == e:
            dPred = 1.0 if pedge._p >= 0.5 else 0.0
            dSeq = 1.0 if e._p >= 30.6 else 0.0
            vals.append((pedge._p, e._p))
            dat.add(pedge._p, dSeq)
            #vals.append( (pedge._p, e._p) )
        else:
            print("CHA!")
            vals.append((0.0, e._p))

    ROCCurve = croc.ROC(dat.sweep_threshold())
    print("ROC area = {0}".format(ROCCurve.area()))
    x,y = zip(*ROCCurve)
    plt.plot(x,y)
    plt.show()
    pred, seq = zip(*vals)
    pred = np.array(pred)
    seq = np.array(seq)
    print("r^2 = {0}".format(np.corrcoef(pred, seq)))
    plt.scatter(pred, seq)
    plt.show()

if __name__ == "__main__":
    main()
