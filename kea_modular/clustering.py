import subprocess,os
from datetime import datetime
def clustering(t,x):
    #Find otu clusters with ANI >99%
    if not os.path.isfile('clusters.csv'):
        print(datetime.now(), 'Clustering initated')
        subprocess.Popen("coverm cluster --genome-fasta-list input_mag_abs_path.tsv --genome-fasta-extension %s --threads %d > clusters.csv" % (x, t), shell=True).wait()
    else:
        print(datetime.now(), 'Cluster file already exists')

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
        if os.path.exists(cluster):
            print(datetime.now(), 'Cluster directory already exists')
        else:
            for key, value in sorted(cluster_dict.items()):
                if len([item for item in value if item]) >= 2: #make sure that there are 2 or more mags in the cluster before running checkm
                    subprocess.Popen("mkdir %s" % cluster, shell=True).wait()
                    genomes = cluster_dict[cluster]
                    for genome in genomes:
                        subprocess.Popen("ln -s %s %s/" % (genome, cluster), shell=True).wait()
                    # running checkm on each cluster
                    print(datetime.now(), 'Quality check initated on  ' + str(cluster))
                    subprocess.Popen("checkm lineage_wf %s %s/checkm --reduced_tree --tab_table -t %d --pplacer_threads %d -x %s -q > %s/checkm.tsv" % (cluster, cluster, t, t, x, cluster), shell=True).wait()
                else:
                    print(datetime.now(), cluster, ' has <2 MAGs, skipping quality check')
                    cluster_dict.pop(key, None) #removes the cluster from the cluster check



    return cluster_dict