from lxml import etree as ET

class CrossValidationSet(object):
	def __init__(self, name):
		self.name = name
		self.edges = []

	def addEdge(self, u, v):
		self.edges.append( sorted((u,v)) )

class CrossValidationTest(object):

	def __init__(self, name):
		self.testSetName = name
		self.cvSets = {}

	def addCVSet(self, cvtest):
		self.cvSets[cvtest.name] = cvtest

def parseCrossValidationTest( fname ):
	T = ET.parse(fname)
	root = T.getroot()

	cvtest = CrossValidationTest( root.attrib["name"] )

	for ts in root.findall("testset"):
		cvset = CrossValidationSet( ts.attrib["name"] )
		for edge in ts.findall("edge"):
			cvset.addEdge( edge.attrib["u"], edge.attrib["v"] )

		cvtest.addCVSet(cvset)

	return cvtest