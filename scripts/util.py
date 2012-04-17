import networkx as nx
import itertools
import collections



def parseHomologyInfo( hfile ):
    '''
    Parse the "Emre" format homology file.  This returns a dictonary
    where the key, k, is an ancestral protein.  The value is a 2-tuple of
    sets.  The first set contains the homologs of k in G1 and the second
    set contains the homologs of k in G2.
    '''
    def parseList( s ) :
        s = s[1:-1]
        return [ str(e)[1:-1] for e in s.split(', ') ] if len(s) > 0 else []

    dupDict = {}
    with open(hfile,'rb') as ifile :
        for l in ifile :
            toks = l.rstrip().split('\t')
            key = toks[0]
            g1Homologs = set([ e for e in parseList(toks[1]) ])
            g2Homologs = set([ e for e in parseList(toks[2]) ])
            dupDict[key] = (g1Homologs, g2Homologs)

    return dupDict

def prepareTree(T, rv, lostNodes):
    '''
    Prepare the tree for various queries we will make later by populating all nodes with
    certain information concerning their subtrees.  In particular, for each node we will
    compute: the leaves beneath it, the set of all subnodes beneath it, and the set of
    extant networks beneath it.
    '''
    # Gather the leaf sets
    for n in T.nodes_iter() :
        T.node[n].update( { 'leaves' : set([]), 'subnodes' : set([]), 'enets' : set([]) } )
    labelWithSubnodes(T, rv, lostNodes)

def labelWithSubnodes( T, root, lostNodes ):
    """
    For each node in the tree T, we add two vertex properties:
    1) leaves -- A set containing all of the leaves rooted at this subtree
    2) subnodes -- A set containing all of the nodes (internal & leaves) in this subtree
    3) enets -- A set containing network identifiers for all extant networks which exist as some leaf in this subtree
    """
    lost = 'LOST'

    def enet(nodes):
        r = set([])
        for n in nodes:
            if n not in lostNodes:
            #   if n.find("LOST") != -1:
            #    r.add(n[:2])
            #else:
                r.add(n[-2:])
            else:
                r.add(lost)
        return r

    for n in nx.dfs_postorder_nodes(T, root) :
        successors = T.successors(n)
        if len(successors) == 0 :
            T.node[n]['subnodes'] =  set([n])
            T.node[n]['leaves'] = T.node[n]['leaves'].union( set([n]) )
            T.node[n]['enets'] = enet([n])
        else :
            subnodes = [ set([n]) ] + [ set([s]) for s in successors ] + [ T.node[s]['subnodes'] for s in successors ]
            for sn in subnodes:
                T.node[n]['subnodes'] = T.node[n]['subnodes'].union( sn )

            leaves = [ T.node[s]['leaves'] for s in successors ]
            T.node[n]['leaves'] = T.node[n]['leaves'].union( *leaves )
            T.node[n]['enets'] = enet( T.node[n]['leaves'] )


def allPairs( nodes ):
    '''
    Yield all pairs of nodes as 2-tuples (including every node with itself)
    '''
    for u,v in itertools.combinations(nodes, 2) :
        if u < v :
            yield u,v
        else :
            yield v, u
    for u in nodes :
        yield u,u

def swap(u,v) :
    '''
    Yup, just swap u and v
    '''
    return v, u

def findRoot( G ) :
    """
    Given a graph G, which is a rooted tree, return the root node.
    """
    import random
    rn = random.choice(G.nodes())
    pred = G.predecessors(rn)
    while len(pred) == 1:
        rn = pred[0]
        pred = G.predecessors(rn)
    return rn

def pairs(lst) :
    n = len(lst)
    return itertools.izip(*[itertools.islice(lst,0,n-1,1), itertools.islice(lst,1,n,1)])

def pathToRoot( n, rv, T ) :
    '''
    Compute the path in the tree T, from node n to the root node, rv.
    '''
    path = []
    if n != rv :
        p = n
        while p != rv :
            path.append(p)
            p = T.predecessors(p)[0]
        path.append(p)
    return path

def edgeCost(edge, direc, cc, dc) :
    '''
    Convenience function to compute the cost of the edge
    between u and v.
    '''
    cost = 0
    u,v = edge

    if dir == 'f' or dir == 'r' :
        cost = cc
    else :
        if u == v :
            cost = cc
        else :
            cost = 2*cc
    return cost

ExistInt = collections.namedtuple('ExistInt', ['creation', 'death'], verbose=False)

def timePenalty( intervalA, intervalB ) :
    if (intervalA.creation <= intervalB.death) and (intervalA.death >= intervalB.creation):
        return 0.0
    else :
        if intervalA.death <= intervalB.creation :
            return intervalB.creation - intervalA.death
        else :
            return intervalA.creation - intervalB.death

def treeToGraph(t) :
    '''
    Convert the Newick format tree, t, into a NetworkX
    DiGraph.  Tree edges will point from parents to children.
    If branch lengths exist in the original tree, they will be
    placed as edge weights in the resulting graph.
    '''
    #nt = t.Children[0]
    #print(nt)
    #t = nt
    #root = nt.Name#t.Name
    root = t.Name
    G = nx.DiGraph()
    # add the root, which exists at time 0 (origin time)
    G.add_node(root)
    G.node[root]['exist interval'] = ExistInt(creation=0, death=0)
    for n in t.preorder():
        if n.Name != root:
            G.add_edge(n.Parent.Name, n.Name, weight=n.Parent.distance(n) )
            # This node comes into existence upon the duplication of it's parent,
            # and exists for the amount of time determined by it's edge length
            parentDeath = G.node[n.Parent.Name]['exist interval'].death
	    etime = float("inf") if n.istip() else n.Parent.distance(n)
            G.node[n.Name]['exist interval'] = ExistInt(creation=parentDeath, death=parentDeath+etime)
    return G
