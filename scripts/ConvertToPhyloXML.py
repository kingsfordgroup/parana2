"""
Usage: ConvertToPhyloXML.py IDIR REGEX

Arguments:
  IDIR     The input directory containing the NXH format trees
  IREGEX   The regular expression (should be in double quotes & python format)
           describing which trees in IDIR should be converted to PhyloXML format.
Options:
  -h       Display this help message.
"""
from docopt import docopt

import re
import os
import subprocess
from progressbar import ProgressBar

def main(options):
	inputDir = options['IDIR']
	regex = options['REGEX']
	inputRE = re.compile(regex)

	print("Converting all trees in {0} matching regex \"{1}\" to PhyloXML format.".format(inputDir, regex))

	cmd="""java -cp ../scripts/forester_1014.jar org.forester.application.phyloxml_converter"""
	for dirname, dirnames, filenames in os.walk(inputDir):
		inputFiles = [ x for x in filenames if inputRE.match(x) ]
		progress = ProgressBar()
		for fn in progress(inputFiles):
			outfn = fn.replace("nwk","xml")
			inputPath = os.path.sep.join([dirname, fn])
			outputPath = os.path.sep.join([dirname, outfn])

			cmdOpts = ["-f=nn", inputPath, outputPath]

			with open(os.devnull,'w') as outstream:
				subprocess.call(cmd.split()+cmdOpts, stdout=outstream, stderr=outstream)

if __name__ == "__main__":
	options = docopt(__doc__, version="ConvertToPhyloXML 1.0")
	main(options)