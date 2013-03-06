# Python modules
import argparse
import os
import logging
from pprint import pprint

import Bio.SeqIO
from Bio.Align.Applications import MuscleCommandline

class MissingDataError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def createParser():
    '''
    This function creates the command line parser.
    '''

    parser = argparse.ArgumentParser(prog="./parana.py")

    parser.add_argument("-g", "--orthoGroups", help="file containing orthology groups")
    parser.add_argument("-d", "--dir", help="directory containing subfolders for genome data")
    parser.add_argument("-o", "--output", help="output directory name")
    return parser


def parseOrthoFile( fname ):
    f = open(fname)
    # first line is a comment
    f.readline()
    species = f.readline().rstrip().split('\t')
    numGroups = int(f.readline().rstrip())
    orthoGroups = {}

    for gn in xrange(numGroups):
        gnum = int(f.readline().rstrip())
        orthoGroups[gnum] = {}
        i = 0
        while len( set(species) - set(orthoGroups[gnum].keys()) ) > 0:
            i += 1
            toks = f.readline().rstrip().split('\t')
            sp = toks[0]
            genes = toks[1:]
            orthoGroups[gnum][sp] = genes
    f.close()

    return (species, orthoGroups)

def getSeqs( seqDict, spec, genes ) :
    import copy
    def getGeneName( s ):
        return s
    def pname( s, genes ):
        return getGeneName(s) in genes

    fgenes = filter( lambda x: pname(x.description, genes), copy.deepcopy(seqDict[spec]) )
    uf = set( [getGeneName(g.description) for g in  fgenes] ) ^ set(genes)
    if len(uf) != 0:
        raise MissingDataError("Could not find genes {0} in species {1}".format(uf,spec))

    for g in fgenes:
	g.description = "{0}_{1}".format(g.description,spec)
	g.id = g.description
    print(fgenes)
    return fgenes

def main():
    '''
    Main method that sets up logging and creates and invokes the command line parser
    '''
    ## Set the logging level
    logging.basicConfig(level=logging.WARNING)

    ## Parse the command line arguments and set the relevant
    ## program options.
    parser = createParser()
    options = parser.parse_args()

    species, orthoGroups = parseOrthoFile( options.orthoGroups )

    specDirs = { sp : os.path.sep.join( [options.dir, sp, "{0}.fasta".format(sp)] ) for sp in species }
    seqDict = {}
    # read in the sequences
    for spec, sd in specDirs.iteritems():
        seqDict[spec] = [seq for seq in Bio.SeqIO.parse( sd, 'fasta' )]

    os.mkdir(options.output)
    ogsDir = os.path.sep.join([options.output, "OrothologyGroupSeqs"])
    os.mkdir( ogsDir )

    for gnum, ogroup in orthoGroups.iteritems():
        seqRecords = []
        for spec, genes in ogroup.iteritems():
            seqRecords += getSeqs( seqDict, spec, genes )
	    print(seqRecords)
        with open( os.path.sep.join([ ogsDir, "{0}.fasta".format(str(gnum)) ]), 'wb' ) as ofile:
            Bio.SeqIO.write( seqRecords, ofile, 'fasta' )

if __name__ == "__main__" : main()
