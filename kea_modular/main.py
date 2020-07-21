#!/usr/bin/env python

import click

from kea_modular.input_processing import input_click
from kea_modular.clustering import clustering
from kea_modular.quality_check import qual_check
from kea_modular.contig_lengths import contig_describe
from kea_modular.nucmer_filtering import filtering


@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', '-o', help='Directory for output polished MAG')
@click.option('-x', default='fa', help='File extension (default: fa)')
@click.option('-t', default=1, help='Number of threads to use (default: 1)')

def full_wf(output, input, x, t):
    input_click(input, output) #Make output directories, process input files
    cluster_dict = clustering(t,x) #Make OTU clusters of input files with ANI >99%
    qual_check(x, cluster_dict) #Quality check using CheckM, removing low quality MAGs, running nucmer
    contig_dict = contig_describe() #Find the contig lengths of the input MAGs to be used for filtering nucmer output
    filtering(contig_dict) #Filter nucmer output to only include end matches


if __name__ == "__main__":
    full_wf()