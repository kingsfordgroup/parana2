"""
Usage: GenerateCrossValidationFile.py --graph=<graph> --folds=<f> --out=<o>

Options
  -h, --help       Display this message
  --graph=<graph>  Input graph
  --folds=<f>      Number of folds to perform (or "loo")
  --out=<o>        Output file
"""
from docopt import docopt
import networkx as nx
import numpy as np
import scipy as sp
from scikits.learn.cross_val import KFold
import random
from lxml import etree as ET
import CrossValidationTools

def main(args):

	G = nx.read_adjlist(args["--graph"])

	leaveOneOut = args["--folds"] == "loo"
	numFolds = None
	if not leaveOneOut:
		numFolds = int(args["--folds"])
	else:
		numFolds = G.size()

	ofname = args["--out"]
	root = ET.Element("cvtest", name="{0}_{1}_test".format(args["--graph"], numFolds))
		
	edges = G.edges()
	random.shuffle(edges)

	kf = KFold(G.size(), numFolds, indices=True)
	for i, (trainIDs, testIDs) in enumerate(kf):
		tset = ET.SubElement(root, "testset", name="fold_{0}".format(i))
		trainEdges = [edges[i] for i in trainIDs]
		testEdges = [edges[j] for j in testIDs]
		for u,v in testEdges:
			ET.SubElement(tset, "edge", u=u, v=v)

	with open(ofname, 'wb') as ofile:
		ofile.write(ET.tostring(root, pretty_print=True))

if __name__ == "__main__":
	arguments = docopt(__doc__, version='GenerateCrossValidationFile v1.0')
	main(arguments)