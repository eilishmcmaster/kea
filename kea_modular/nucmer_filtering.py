import os

def filtering(cluster_dict,  sorted_dict):
    for cluster in cluster_dict:
        os.chdir(cluster)
        #repmag is first in sorted_dict

        #find the repmag for that cluster and find its contig ends
        #keep matches with 1-200bp or length-200
        for mag in sorted_dict:
            with open(nucmer.coords, ) as nucmer_file:
                #find dictionary entry corresoponding to alignment

        os.chdir('../')