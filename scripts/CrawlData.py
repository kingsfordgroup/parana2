import urllib2
import networkx as nx
import sys

def firstEnt(s):
  loc = s.find(">", 2)
  return s[:loc]

def main():
  gname = sys.argv[1]  
  species = sys.argv[2]
  ofname = sys.argv[3]
  ofile = open(ofname,'wb')

  G = nx.read_adjlist(gname)
  notFound = []
  fs = 'http://www.uniprot.org/uniprot/?query={0}+AND+organism%3A{1}&sort=score&format=fasta&limit=3' 
  for i,n in enumerate(G.nodes()):
    try:
      on = float(n)
      on = "ORF"+n  
    except Exception as e:
      on = n

    url = fs.format(on, species)
    print("fetching {0} using {1}".format(on,url)) 
    req = urllib2.urlopen( url ) 
    e = firstEnt( req.read() ) 
    if e == "":
      notFound.append(on)
      #raise NameError("Could not find {0} @ {1}".format(n, url))
    else:
      ofile.write( ">{0}\n".format(n)+"\n".join(e.split("\n")[1:]) )

  ofile.close()
  print("couldn't find {0}".format(notFound))

if __name__ == "__main__": main()

