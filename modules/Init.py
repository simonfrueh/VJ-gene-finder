import os
import argparse

# parse command line arguments
def parse_arguments():
	
	# Initialize command line argument parser
	parser = argparse.ArgumentParser(
		prog = "VJ-gene-finder",
		description = "Extraction of chicken TCR VJ gene segments" + 
					  "from a reference genome",
		epilog = "https://github.com/simonfrueh/VJ-gene-finder")
	 
	# Add and parse arguments
	# To Do: option for Output directory
	parser.add_argument("filename", help = "Input file")
	args = parser.parse_args()

	if args.filename:
		print("Input file: %s" % args.filename)

	# parse path
	path = os.path.dirname(args.filename).replace(" ","")
	
	# Add trailing slahs for later usage
	if len(path) > 0: 
		if path[-1] != "/": 
			path += "/"
	
	# parse filename
	filename = os.path.basename(args.filename) 
	
	return path, filename
