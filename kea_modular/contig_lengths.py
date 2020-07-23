import subprocess
import pandas as pd

def contig_describe():
    # run the perl describe contig script on each input MAG and send output to CSV
    with  open('input_mag_abs_path.tsv', 'r') as input_mags:
        with open('contig_describe.csv', 'w') as csvfile:
            for mag in input_mags:
                perl = subprocess.Popen(["perl", "../../kea_modular/describe_fasta.pl", mag], stdout=subprocess.PIPE)
                perl_output = perl.communicate()[0].decode('UTF-8')
                for line in perl_output:
                    csvfile.write("%s" % line)

        # read in CSV we just made
        col_list = ['Number', 'Name', 'Length']
        perl_csv = pd.read_csv('contig_describe.csv', sep='\t', usecols=col_list)
        perl_csv = perl_csv[perl_csv.Number != 'Number']
        perl_csv = perl_csv.drop(['Number'], axis=1)

        # makes dictionary that looks like {'SHRC1_26_bin_22_contig_968': ['8341'], 'SHRC1_26_bin_22_contig_9697': ['14596']}
        contig_dict = perl_csv.applymap(str).groupby('Name')['Length'].apply(str).to_dict()

    return contig_dict



