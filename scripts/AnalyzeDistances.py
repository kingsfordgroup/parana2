import numpy as np

def getExistenceIntervals(t):

    eints = {}
    rname = t.root().Name
    eints[ rname ] = (-np.inf, 0)

    for n in t.preorder():
        name = n.Name
        if name != rname:
            eints[name] = ( eints[ n.Parent.Name ][1], eints[ n.Parent.Name ][1] + n.Parent.distance( n )  )

    return eints

def dist( n1, n2, eints ) :
    if ( eints[n1][1] >= eints[n2][0] and eints[n2][1] >= eints[n1][0] ):
        return 0.0
    elif ( eints[n1][0] > eints[n2][1] ):
        return eints[n1][0] - eints[n2][1]
    elif ( eints[n2][0] > eints[n1][1] ):
        return eints[n2][0] - eints[n1][1]
    else:
        return np.inf

def main():
    pass

if __name__ == "__main__":
    main()