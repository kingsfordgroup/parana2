"""
Usage: ClustalToNexus.py --in=<in> --out=<out> [--infmt=<infmt>] [--outfmt=<outfmt>]

Options:
  --in=<in>     Input file
  --out=<out>   Output file
  --infmt=<i>   Input format [default: clustal]
  --outfmt=<o>  Input format [default: nexus]
"""
from docopt import docopt
import Bio
from Bio import SeqIO

def main(args):
	ifname = args["--in"]
	ofname = args["--out"]
	infmt = args['--infmt']
	outfmt = args['--outfmt']

	seqs = []
	alpha = Bio.Alphabet.ProteinAlphabet()
	with open(ifname) as ifile:
		seqs = [s for s in SeqIO.parse(ifile,infmt,alpha)]

	with open(ofname,'wb') as ofile:
		SeqIO.write(seqs, ofile, outfmt)

if __name__ == "__main__":
  arguments = docopt(__doc__, version='ClustalToNexus v1.0')
  main(arguments)
