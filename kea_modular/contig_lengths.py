import subprocess
import pandas as pd

def contig_describe():
    # run the perl describe contig script on each input MAG and send output to CSV
    with  open('input_mag_abs_path.tsv', 'r') as input_mags:
        with open('../contig_describe.csv', 'w') as csvfile:
            for mag in input_mags:
                perl = subprocess.Popen(["perl", "../../kea_modular/describe_fasta.pl", mag], stdout=subprocess.PIPE)
                perl_output = perl.communicate()[0].decode('UTF-8')
                for line in perl_output:
                    if not line.startswith("Number"):
                        csvfile.write("%s" % line)

    # read in CSV we just made
    col_list = ['Number', 'Name', 'Length']
    perl_csv = pd.read_csv('../contig_describe.csv', sep='\t')
    #perl_csv = perl_csv[perl_csv.Number != 'Number']
    perl_csv = perl_csv.drop(['Number'], axis=1)
    # makes dictionary that looks like {'NODE_111_length_94059_cov_6.576444': {'Length': '94059'}, 'NODE_131_length_88844_cov_5.932694': {'Length': '88844'}}
    contig_dict = perl_csv.set_index('Name').T.to_dict()
    # refer like this: print(contig_dict['NODE_131_length_88844_cov_5.932694']['Length'])
    return contig_dict
