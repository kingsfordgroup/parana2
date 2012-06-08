import logging
import networkx as nx
from util import allPairs

def calcStats(GX,GA) :
    '''
    Calculate precision, recall and f1 statistics. GX is the "true" ancestor network,
    and GA is our prediction.
    '''
    tp = set( [ e for e in GX.edges() if GA.has_edge(e[0],e[1]) ] )
    fp = set( [ e for e in GA.edges() if not GX.has_edge(e[0], e[1]) ] )
    fn = set( [ e for e in GX.edges() if not GA.has_edge(e[0], e[1]) ] )

    precision = float( len(tp) ) / (len(tp)+len(fp)) if (len(tp)+len(fp)) > 0 else 0.0
    recall =  float(len(tp)) / (len(tp) + len(fn)) if (len(tp)+len(fn)) > 0 else 0.0
    prsum = (precision + recall)
    f1 = 2.0 * (precision * recall) / (precision + recall) if prsum > 0 else 0.0

    logging.info("===== False Pos =====")
    logging.info(fp)
    logging.info("===== False Neg =====")
    logging.info(fn)

    return precision, recall, f1

def printStats( recon, orig ) :
    '''
    Given a reconstructed ancestral network (recon) and a __ground truth__ network (orig),
    compute the precision, recall and F1-score of the reconstruction.
    '''
    tp = [(u,v) for u,v in recon.edges_iter() if orig.has_edge(u,v)]
    tn = [(u,v) for u,v in allPairs( orig.nodes() ) if not (recon.has_edge(u,v) and orig.has_edge(u,v)) ]
    fp = [(u,v) for u,v in recon.edges_iter() if not orig.has_edge(u,v)]
    fn = [(u,v) for u,v in orig.edges_iter() if not recon.has_edge(u,v)]
    prec = len(tp) / float( len(tp) + len(fp) ) if float( len(tp) + len(fp) ) > 0 else 0.0
    sens = rec = len(tp) / float( len(tp) + len(fn) ) if float( len(tp) + len(fp) ) > 0 else 0.0
    print("Precision = {0}, Recall = {1}, F1-Score = {2}".format( prec, rec, 2*(prec*rec)/(prec+rec) if (prec+rec) > 0 else 0.0 ) )
    omspec = 1 - ( len(tn) / float(len(tn)+len(fp)) )
    return (omspec, sens)
