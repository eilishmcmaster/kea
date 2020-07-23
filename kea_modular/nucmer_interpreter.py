from glob import iglob
import pandas as pd
import os

def nucmer_interpreter(rep_mag, mag):
    #change coord file to csv
    coordfile = ('%s.coords' % (rep_mag + '_vs_' + mag))
    base = os.path.splitext(coordfile)[0]
    os.rename(coordfile, base + '.csv')

    colnames=['ref_start','mag_start','length_1','%_identity','ref_contig']
    pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.width', None, 'display.max_colwidth', -1)
    data = pd.read_csv(next(iglob('%s.csv' % (rep_mag + '_vs_' + mag))), sep='|', header=None, names=colnames)

    #remove the top junk lines
    data.drop([0,1,2,3], axis=0, inplace=True)

    #split the double columns (everything is converted to strings)
    data[['ref_contig','other_contig']]=data['ref_contig'].str.split('\t', expand=True)
    data[['ref_start','ref_end']]=data['ref_start'].str.split(expand=True)
    data[['mag_start','mag_end']]=data['mag_start'].str.split(expand=True)
    data[['length_1','length_2']]=data['length_1'].str.split(expand=True)

    return data
