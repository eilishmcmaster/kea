import pandas as pd

def gap_overlap(nucmer_array, contig_dict):
    for index, alignment in nucmer_array.iterrows():
        ref_contig = alignment['ref_contig']
        ref_contig = str(ref_contig[1:])
        ref_contig_length = contig_dict[ref_contig]['Length']
        ref_contig_length = int(ref_contig_length)

        other_contig = alignment['other_contig']
        other_contig_length = contig_dict[other_contig]['Length']
        other_contig_length = int(other_contig_length)

        if int(alignment['mag_start']) <=200 and int(alignment['ref_end']) >= ref_contig_length-200 or\
                int(alignment['ref_start']) <= 200 and int(alignment['mag_end']) >= other_contig_length -200:
            print('Alignment removed')
            nucmer_array.remove(alignment)
        else:
            print('Alignment approved')
    return nucmer_array


