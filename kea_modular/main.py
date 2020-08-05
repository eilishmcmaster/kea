#!/usr/bin/env python

import click

from input_processing import input_click
from clustering import clustering
from nucmer import nucmer
from contig_lengths import contig_describe


@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', '-o', help='Directory for output polished MAG')
@click.option('-x', default='fa', help='File extension (default: fa)')
@click.option('-t', default=1, help='Number of threads to use (default: 1)')

def full_wf(output, input, x, t):
    input_click(input, output) #Make output directories, process input files
    cluster_dict = clustering(t,x) #Make OTU clusters of input files with ANI >99%
    contig_dict = contig_describe() #Find the contig lengths of the input MAGs to be used for filtering nucmer output
    nucmer(x, cluster_dict, contig_dict) #Quality check using CheckM, removing low quality MAGs, running nucmer

if __name__ == "__main__":
    full_wf()