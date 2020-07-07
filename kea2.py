#!/usr/bin/env python

import os
import subprocess
import click
from os import path

#python kea2.py -i /mnt/c/Users/Eilish\ McMaster/PycharmProjects/test_files -o output_dir -x
#python 02072020.py -i ~/kea_wd/test_files/ -o out_dir -x fa -t 20


@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')
@click.option('--output', '-o', help='Directory for output polished MAG')
@click.option('-x', default='fa', help='File extension (default: fa)')
@click.option('-t', default=1, help='Number of threads to use (default: 1)')
def full_wf(output, input, x, t):

    #Make output directory as specified by --output
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

    #Make a list of all of the input files with their absolute path to use later
    files_in_input_dir = os.listdir(input)
    input_mag_abs_path = open('start.tsv', 'w')
    for input_file in files_in_input_dir:
        input_mag_abs_path.write("%s/%s\n" % (input, input_file))
    input_mag_abs_path.close()
    #remove spaces
    fin = open('start.tsv', 'rt')
    fout = open('input_mag_abs_path.tsv', 'wt')
    for line in fin:
        fout.write(line.replace(' ', '\ '))
    fin.close()
    fout.close()
    print('File containing list of input files with absolute path created (input_mag_abs_path.tsv)')

    #Find otu clusters with ANI >99%
    print('Clustering initated')
    subprocess.Popen("coverm cluster --genome-fasta-list input_mag_abs_path.tsv --genome-fasta-extension %s --threads %d > clusters.csv" % (x, t), shell=True).wait()

    #Make dictionary of different clusters
    from collections import defaultdict
    cluster_dict = defaultdict(set)
    name_dict = {}
    with open('clusters.csv', 'r+') as data:
        # Cluster identifier
        cluster_id = 1
        for line in data:
            key, val = line.strip().split()
            # Add to name dict
            if key not in name_dict.keys():
                name_dict[key] = "cluster_" + str(cluster_id)
                cluster_id += 1
            # Gets cluster name from name_dict, appends new cluster name to cluster_dict if it does not already exist
            # and initiates new cluster set
            # If it does exist, appends to new value to cluster set
            cluster_dict[name_dict[key]].add(val)


    # #checkm on all input files
    # print('Quality check initated')
    # subprocess.Popen("checkm lineage_wf %s checkm_output --reduced_tree --tab_table -t %d --pplacer_threads %d -x %s -q >checkm.tsv" % (input, t, t, x), shell=True).wait()

    #Run checkM on individual clusters


if __name__ == "__main__":
    full_wf()


