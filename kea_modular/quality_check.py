import pandas as pd
import json, subprocess, os
import numpy as np

def qual_check(x, cluster_dict):
    for cluster in cluster_dict:
        #make dictionary for checkm.tsv
        col_list = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4','5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        checkm_file = pd.read_csv("%s/checkm.tsv" % cluster, sep='\t', usecols=col_list)

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

        #remove bad quality MAGs from clusters
        for index, row in checkm_file.iterrows():
            if row['Contamination'] >= 5:
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                del combined_checkm[row['Bin Id']]
                print(str(row['Bin Id']) + ' removed due to high contamination (>=5%)')
            if row['Completeness'] > checkm_file['Completeness'].mean() + np.std(checkm_file['Completeness']) or row['Completeness'] < checkm_file['Completeness'].mean() - np.std(checkm_file['Completeness']):
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                del combined_checkm[row['Bin Id']]
                print(row['Bin Id'] + ' removed due to uncharacteristic completeness (contamination outside 1 standard deviation of the cluster mean)')

        full_checkm_dict = {}
        # calculate genome score for each MAG within the cluster
        for mag, items in combined_checkm.items():
            genome_score = items['Completeness'] - 5 * items['Contamination'] - 5 * (items['# contigs'] / 100) - 5 * (items['# ambiguous bases'] / 10000)
            combined_checkm[mag]['Genome score'] = genome_score
            print(mag, ' genome score: ', items['Genome score'])
            full_checkm_dict[mag] = items['Genome score']

        # sort MAGs by genome score
        import operator
        sorted_dict = dict(sorted(full_checkm_dict.items(), key=operator.itemgetter(1), reverse=True))

        # find the highest MAG (individually) and remove it from the dictionary
        rep_mag = max(full_checkm_dict.items(), key=operator.itemgetter(1))[0]
        print(rep_mag, ' has the highest genome score for cluster')
        print(sorted_dict)

        os.chdir(cluster)
        # use nucmer to compare representative mag to other mags in genome score order
        for mag in sorted_dict:
            subprocess.Popen("nucmer --prefix=%s %s %s --coords" % (rep_mag + '_vs_' + mag, rep_mag + '.' + x, mag + '.' + x), shell=True).wait()
            print(rep_mag + ' aligned with ' + mag + ' in ' + cluster)
        subprocess.Popen("rm %s*" % (rep_mag + '_vs_' + mag), shell = True).wait()
        os.chdir('../')

    return sorted_dict
