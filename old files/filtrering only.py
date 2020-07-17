#!/usr/bin/env python

import os
import subprocess
import click
import pandas as pd

os.chdir('out_dir/kea_wd/')
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

# making directories for each cluster
for cluster in cluster_dict:
    col_list = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4',
                '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
    file = pd.read_csv("%s/checkm.tsv" % (cluster), sep='\t', usecols=col_list)

    for index, row in file.iterrows():
        if row['Contamination'] >= 5:
            subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
            print(str(row['Bin Id']) + 'removed due to high contamination (>=5%)')