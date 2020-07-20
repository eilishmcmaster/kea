#!/usr/bin/env python

import click

from kea_modular.input_processing import input_click
from kea_modular.clustering import clustering
from kea_modular.quality_check import qual_check
from kea_modular.nucmer import sort_n_nucmer

@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', '-o', help='Directory for output polished MAG')
@click.option('-x', default='fa', help='File extension (default: fa)')
@click.option('-t', default=1, help='Number of threads to use (default: 1)')

def full_wf(output, input, x, t):
    input_click(input, output) #Make output directories, process input files
    cluster_dict = clustering(t,x) #Make OTU clusters of input files with ANI >99%
    combined_checkm = qual_check(x, cluster_dict) #Quality check using CheckM and removing low quality MAGs
    sort_n_nucmer(cluster_dict, combined_checkm,x) #Sorting MAGs by genome score, designating representative MAG and running nucmer


if __name__ == "__main__":
    full_wf()