#!/usr/bin/env python

import os
import subprocess
import click
from os import path

#python using_coverm.py --inputs /mnt/c/Users/Eilish\ McMaster/PycharmProjects/kea/test_files --output output_dir cluster
#python using_fastani.py --inputs ~/kea_wd/test_files --output output_dir cluster


@click.group()
@click.version_option()
@click.option('--inputs', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', help='Directory for output polished MAG')
#@click.option('--output', type=click.Path(), help='Directory for output polished MAG')
#@click.option('--extension', default='fa', help='File extension (fa, fna, fasta)')
#@click.option('--threads', default=1, help='Number of threads to use')
def setup(output, inputs):
    '''Sets up the directories and processes input'''
    #make output directory as specified by --output
    if os.path.exists(output):
        print('Error: Output directory already exists')
        exit()
    else:
        os.makedirs(output)
        os.chdir(output)
        print('Output directory created')
        #make a working directory
        os.makedirs('kea_wd')
        os.chdir('kea_wd')
        print('Working directory created')


    #make a list of all of the input files with their absolute path to use later
    files_in_input_dir = os.listdir(inputs)
    input_mag_abs_path = open('input_mag_abs_path.tsv', 'w')
    for input_file in files_in_input_dir:
        input_mag_abs_path.write("%s/%s\n" % (inputs, input_file))
    input_mag_abs_path.close()
    print('File containing list of input files with absolute path created')

@setup.group()
def bullshit():
    '''fuck you'''
@bullshit.command()
@click.option('--extension', default='fa', help='File extension (fa, fna, fasta)')
@click.option('--threads', default=1, help='Number of threads to use')
def cluster(extension, threads):
    '''Uses coverm to determine >99% ANI OTU clusters'''
    click.echo(threads)
    click.echo(extension)
    #subprocess.Popen("coverm cluster --genome-fasta-list input_mag_abs_path.tsv --genome-fasta-extension %s --threads %d > clusters.tsv", (shell=True)).wait()



if __name__ == "__main__":
    setup()
    cluster()
    #coverm_clustering()
