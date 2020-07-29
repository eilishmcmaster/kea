from glob import iglob
import pandas as pd
import os

def nucmer_interpreter(rep_mag, mag):
    #change coord file to csv
    coordfile = ('%s.coords' % (rep_mag + '_vs_' + mag))
    base = os.path.splitext(coordfile)[0]
    os.rename(coordfile, base + '.csv')

    colnames=['ref_start','mag_start','length_1','identity','ref_contig']
    pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.width', None, 'display.max_colwidth', None)
    data = pd.read_csv(next(iglob('%s.csv' % (rep_mag + '_vs_' + mag))), sep='|', header=None, names=colnames)

    #remove the top junk lines
    data.drop([0,1,2,3], axis=0, inplace=True)

    #split the double columns (everything is converted to strings)
    data[['ref_contig','other_contig']]=data['ref_contig'].str.split('\t', expand=True)
    data[['ref_start','ref_end']]=data['ref_start'].str.split(expand=True)
    data[['mag_start','mag_end']]=data['mag_start'].str.split(expand=True)
    data[['length_1','length_2']]=data['length_1'].str.split(expand=True)

    # convert appropriate columns to integers
    data["length_1"] = data["length_1"].astype(int)
    data["length_2"] = data["length_2"].astype(int)
    data["ref_start"] = data["ref_start"].astype(int)
    data["mag_start"] = data["mag_start"].astype(int)
    data["ref_end"] = data["ref_end"].astype(int)
    data["mag_end"] = data["mag_end"].astype(int)
    data["identity"] = data["identity"].astype(float)

    # remove less good matches
    data = data[data.length_1 > 100]
    data = data[data.length_2 > 100]
    data = data[data.identity > 97]
    nucmer_array=data
    return nucmer_array
