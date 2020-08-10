import pandas as pd
import json, subprocess, os
import numpy as np
from nucmer_interpreter import nucmer_interpreter
from selecting_gap_overlaps import gap_overlap
from remove_duplicates import remove_dupes
from align_and_consensus import consensus_maker
from datetime import datetime

def nucmer(x, cluster_dict, contig_dict, assigned_mag):
    for cluster in cluster_dict:
        print(datetime.now(), ' Making dictionary from checkM output for '+ cluster)
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
        print(datetime.now(), ' Removing bad quality MAGs from ' + cluster)
        for index, row in checkm_file.iterrows():
            if row['Contamination'] >= 5:
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                del combined_checkm[row['Bin Id']]
                print(datetime.now(), (row['Bin Id']) + ' removed due to high contamination (>=5%)')
            if row['Completeness'] > checkm_file['Completeness'].mean() + np.std(checkm_file['Completeness']) or row['Completeness'] < checkm_file['Completeness'].mean() - np.std(checkm_file['Completeness']):
                subprocess.Popen("rm %s/%s.%s" % (cluster, row['Bin Id'], x), shell=True).wait()
                del combined_checkm[row['Bin Id']]
                print(datetime.now(), row['Bin Id'], ' removed due to uncharacteristic completeness (contamination outside 1 standard deviation of the cluster mean)')

        full_checkm_dict = {}
        # calculate genome score for each MAG within the cluster
        print(datetime.now(), ' Calculating genome score for MAGs in ' + cluster)
        for mag, items in combined_checkm.items():
            genome_score = items['Completeness'] - 5 * items['Contamination'] - 5 * (items['# contigs'] / 100) - 5 * (items['# ambiguous bases'] / 10000)
            combined_checkm[mag]['Genome score'] = genome_score
            full_checkm_dict[mag] = items['Genome score']

        # sort MAGs by genome score
        import operator
        sorted_dict = dict(sorted(full_checkm_dict.items(), key=operator.itemgetter(1), reverse=True))

        #determining if assigned mag is in cluster
        in_cluster = False
        for key, item in sorted_dict.items():
            if key == assigned_mag:
                in_cluster = True
                print(datetime.now(), 'Assigned MAG is in cluster ', cluster)
        # find the highest MAG (individually) and remove it from the dictionary
        if assigned_mag == '' or in_cluster == False:
            print(datetime.now(), 'Representative MAG not assigned, or not found in ', cluster)
            rep_mag = max(full_checkm_dict.items(), key=operator.itemgetter(1))[0]
            print(datetime.now(), rep_mag + ' has the highest genome score for cluster')
        else:
            print(datetime.now(), assigned_mag, ' designated as representative MAG')
            rep_mag = assigned_mag

        if rep_mag in sorted_dict:
            del sorted_dict[rep_mag]

        os.chdir(cluster)

        # use nucmer to compare representative mag to other mags in genome score order
        for mag in sorted_dict:
            subprocess.Popen("nucmer --prefix=%s %s %s --coords -q" % (rep_mag + '_vs_' + mag, rep_mag + '.' + x, mag + '.' + x), shell=True).wait()
            print(datetime.now(), rep_mag + ' aligned with ' + mag + ' in ' + cluster)
            subprocess.Popen('rm *delta', shell=True).wait()

            #nucmer array is the polished nucmer output with matches >100bp and identity >97% of rep mag vs other mag
            print(datetime.now(), 'Interpreting nucmer output')
            nucmer_array = nucmer_interpreter(rep_mag, mag)
            if not nucmer_array.empty:
                print(datetime.now(), 'Removing non-end alignments')
                filtered_nucmer = gap_overlap(nucmer_array,contig_dict)
                if not filtered_nucmer.empty:
                    print(datetime.now(), 'Removing duplicates >2')
                    final_nucmer = remove_dupes(filtered_nucmer)
                    if not final_nucmer.empty:
                        print(datetime.now(), 'Closing ', rep_mag, ' gaps using ', mag)
                        results = consensus_maker(final_nucmer, mag, rep_mag,x, contig_dict)
                        rep_mag = results[0]
                        contig_dict = results[1]
                        print(datetime.now(), 'Improvement complete for this pair')
            else:
                print(datetime.now(), 'Array empty after filtering, moving on to next MAG')


        os.chdir('../')