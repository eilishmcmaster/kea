# Make dictionary of different clusters
from collections import defaultdict
import numpy as np
import pandas as pd
import subprocess
import json

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

# path = 'output_dir/cluster_11/checkm/storage/bin_stats.analyze.tsv'
#
# with open(path) as f:
#     data=json.load(f)
# print(data)



#works
col_list = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
file = pd.read_csv("checkm.tsv", sep='\t', usecols=col_list)
print(file['Contamination'])

for cluster in cluster_dict:
    for index, row in file.iterrows():
        # if row['Contamination'] >= 5:
        #     subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], fa), shell=True).wait()
        #     print(row['Bin Id'])
        if row['Completeness'] > file['Completeness'].mean() + np.std(file['Completeness']) or row['Completeness'] < file['Completeness'].mean() - np.std(file['Completeness']):
            print(row['Bin Id'])

# stdev = ('Standard deviation of contamination: ' + str(np.std(file['Contamination'])) + '\n')
# print(stdev)

# f = open('checkm.tsv', "r")
# lines = f.readlines()
# print(lines)
# result = []
# for x in lines:
#     result.append(x.split(' ')[1])
# f.close()
# print(result)

# for cluster in cluster_dict:
#     file = open('checkm.tsv', 'a+')
#     file.write('Stdev' + str(np.std(numbers)) + "\n")
