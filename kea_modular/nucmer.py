import subprocess, os

def sort_n_nucmer(cluster_dict, combined_checkm, x):
    for cluster in cluster_dict:

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
        del sorted_dict[rep_mag]

        os.chdir(cluster)
        # use nucmer to compare representative mag to other mags in genome score order
        for mag in sorted_dict:
            subprocess.Popen("nucmer --prefix=%s %s %s --coords" % (rep_mag + '_vs_' + mag, rep_mag + '.' + x, mag + '.' + x), shell=True).wait()
            print(rep_mag + ' aligned with' + mag + ' in ' + cluster)

        os.chdir('../')