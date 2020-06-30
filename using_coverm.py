#!/usr/bin/env python

import os
import subprocess
import click

#python using_fastani.py --inputs /mnt/c/Users/Eilish\ McMaster/PycharmProjects/kea/test_files --output output_dir cluster
#python using_fastani.py --inputs ~/kea_wd/test_files --output output_dir cluster


@click.group()
@click.option('--inputs', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', help='Directory for output polished MAG')
#@click.option('--output', type=click.Path(), help='Directory for output polished MAG')
#@click.option('--extension', default='fa', help='File extension (fa, fna, fasta)')
#@click.option('--threads', default=1, help='Number of threads to use')

def setup(output, inputs):
    #make output directory as specified by --output
    os.makedirs(output)
    os.chdir(output)
    #make a working directory
    os.makedirs('kea_wd')
    os.chdir('kea_wd')

    #make a list of all of the input files with their absolute path to use later
    files_in_input_dir = os.listdir(inputs)
    input_mag_abs_path = open('input_mag_abs_path.tsv', 'w')
    for input_file in files_in_input_dir:
        input_mag_abs_path.write("%s/%s\n" % (inputs, input_file))
    input_mag_abs_path.close()


@setup.command()
def cluster():
    subprocess.Popen("coverm cluster --genome-fasta-list input_mag_abs_path.tsv --genome-fasta-extension fa --threads 20 > clusters.tsv" , shell=True).wait()

if __name__ == "__main__":
    setup()
    cluster()
