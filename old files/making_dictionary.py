# #making cluster dictionary
import csv
import os
import subprocess
from collections import defaultdict


cluster_dict = defaultdict(set)
name_dict = {}
with open('clusters.csv', 'r+') as data:
    # Cluster identifier
    cluster_id = 1
    for line in data:
        key, val = line.strip().split("\t")
        # Add to name dict
        if key not in name_dict.keys():
            name_dict[key] = "cluster_" + str(cluster_id)
            cluster_id += 1
        # Gets cluster name from name_dict, appends new cluster name to cluster_dict if it does not already exist
        # and initiates new cluster set
        # If it does exist, appends to new value to cluster set
        cluster_dict[name_dict[key]].add(val)
print(cluster_dict)

# making directories for each cluster
for cluster in cluster_dict:
    subprocess.Popen("mkdir %s" % (cluster), shell=True).wait()
    genomes = cluster_dict[cluster]
    for genome in genomes:
        genome = genome.replace("/mnt/c/Users/Eilish\ McMaster/PycharmProjects/", "../../")
        print(genome)
        subprocess.Popen("ln -s %s %s/" % (genome, cluster), shell=True).wait()

    # os.chdir(cluster)
    # subprocess.Popen("ln -s %s" % (cluster_dict(key)), shell=True).wait()


# from collections import defaultdict
# clusters = open('dict_practice.tsv', 'r')
# cluster_dict = defaultdict(list)
# for line in clusters:
#     cluster_dict[line[0:]].append(line[1:])
#
# print(cluster_dict)
#
# cluster_dict = {}
#
# with open('clusters.csv', 'r') as file:
#     lst = []
#     for line in file:
#         lst.append([str(x) for x in line.split('\t')])
#
#
# column1=[x[0] for x in lst ]
# column2=[x[1] for x in lst]
#
# #print(lst)
# print(column1)
# print(column2)













#this works
# d={}
# with open('clusters.csv', 'r') as f:
#     for line in f:
#         (key, val) = line.split('\t')
#         d[str(key)] =val
#
# print(d.get("/srv/home/s4477970/kea_wd/test_files/2012_S3D.11.fa"))


