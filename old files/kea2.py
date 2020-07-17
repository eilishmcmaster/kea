#!/usr/bin/env python

import os
import subprocess
import click
import pandas as pd
from collections import defaultdict
import numpy as np
from os import path
import json

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
            key, val = line.strip().split('\t')
            # Add to name dict
            if key not in name_dict.keys():
                name_dict[key] = "cluster_" + str(cluster_id)
                cluster_id += 1
            # Gets cluster name from name_dict, appends new cluster name to cluster_dict if it does not already exist
            # and initiates new cluster set
            # If it does exist, appends to new value to cluster set
            cluster_dict[name_dict[key]].add(val)

    # making directories for each cluster
    for cluster in cluster_dict:
        subprocess.Popen("mkdir %s" % (cluster), shell=True).wait()
        genomes = cluster_dict[cluster]
        for genome in genomes:
            #genome = genome.replace("/mnt/c/Users/Eilish\ McMaster/PycharmProjects/", "../../")
            print(genome)
            subprocess.Popen("ln -s %s %s/" % (genome, cluster), shell=True).wait()
        # running checkm on each cluster
        print('Quality check initated on  ' + str(cluster))
        subprocess.Popen("checkm lineage_wf %s %s/checkm --reduced_tree --tab_table -t %d --pplacer_threads %d -x %s -q > %s/checkm.tsv" % (cluster, cluster, t, t, x, cluster), shell=True).wait()

        #make dictionary for checkm.tsv
        col_list = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4',
                    '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        checkm_file = pd.read_csv("%s/checkm.tsv" % cluster, sep='\t', usecols=col_list)

        #remove bad quality MAGs from clusters
        for index, row in checkm_file.iterrows():
            if row['Contamination'] >= 5:
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                print(str(row['Bin Id']) + ' removed due to high contamination (>=5%)')
            if row['Completeness'] > checkm_file['Completeness'].mean() + np.std(checkm_file['Completeness']) or row['Completeness'] < checkm_file['Completeness'].mean() - np.std(checkm_file['Completeness']):
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                print(row['Bin Id'] + ' removed due to uncharacteristic completeness (contamination outside 1 standard deviation of the cluster mean)')

        #turn checkm_file into dictionary
        checkm_dict = checkm_file.to_dict(orient="records")
        checkm_tsv = {}
        for mag in checkm_dict:
            checkm_tsv[mag["Bin Id"]] = mag

        #make dictionary for Checkm QA
        checkm_qa = {}
        try:
            with open("%s/checkm/storage/bin_stats.analyze.tsv" % cluster) as qa:
                for line in qa:
                    mag, items = line.strip().split("\t")
                    items = items.replace("'", '"')
                    items = json.loads(items)

                    checkm_qa[mag] = items
        except FileNotFoundError:
            print("%s stats not found" % cluster)

        from deepmerge import always_merger

        #make combined checkm and checkm QA dictionary
        combined_checkm = always_merger.merge(checkm_tsv, checkm_qa)

        stupid_dict = {}
        # calculate genome score for each MAG within the cluster
        for mag, items in combined_checkm.items():
            genome_score = items['Completeness'] - 5 * items['Contamination'] - 5 * (items['# contigs'] / 100) - 5 * (items['# ambiguous bases'] / 10000)
            combined_checkm[mag]['Genome score'] = genome_score
            print(mag, ' genome score: ', items['Genome score'])

            stupid_dict[mag] = items['Genome score']

        # sort MAGs by genome score
        print(stupid_dict)
        import operator
        sorted_dict = dict(sorted(stupid_dict.items(), key=operator.itemgetter(1), reverse=True))
        print(sorted_dict)

        # find the highest MAG (individually) and remove it from the dictionary
        rep_mag = max(stupid_dict.items(), key=operator.itemgetter(1))[0]
        print(rep_mag, ' has the highest genome score for cluster')
        del sorted_dict[rep_mag]
        print(sorted_dict)

        #use nucmer to compare representative mag to other mags in genome score order
        for mag in sorted_dict:
            # subprocess.Popen("nucmer --prefix=%s %s/%s %s/%s --coords" % (rep_mag + '_vs_' + mag, cluster, rep_mag, cluster, mag] ), shell=True).wait()
            print('rep_vs '+ cluster, rep_mag, cluster, mag)

    


if __name__ == "__main__":
    full_wf()


