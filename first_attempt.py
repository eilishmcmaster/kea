#!/usr/bin/python

# this is my first go at doing stuff


import argparse


def run(args):
	filename = args.input # these match the "dest": dest="input"
	output_filename = args.output # from dest="output"
	qual = args.quality_score # default is I

	# Do stuff


def main():
	parser=argparse.ArgumentParser(description="Convert a fastA file to a FastQ file")
	parser.add_argument("-in",help="fasta input file" ,dest="input", type=str, required=True)
	parser.add_argument("-out",help="fastq output filename" ,dest="output", type=str, required=True)
	parser.add_argument("-q",help="Quality score to fill in (since fasta doesn't have quality scores but fastq needs them. Default=I" ,dest="quality_score", type=str, default="I")
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()