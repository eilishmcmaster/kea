#!/usr/bin/env python
import os
import subprocess
import click

#python using_fastani.py --inputs /mnt/c/Users/Eilish\ McMaster/PycharmProjects/kea/test_files --output output_dir

@click.command()
@click.option('--output', type=click.Path(), help='Directory for output polished MAG')
@click.option('--inputs', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')

#@click.option('--threads', default=16, help='Number of threads to use')

#make output directory
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
        input_mag_abs_path.write("%s \n" % os.path.abspath(input_file))
    input_mag_abs_path.close()


if __name__ == "__main__":
    setup()