import pandas as pd

def gap_overlap(nucmer_array, contig_dict):
    col_names = ['ref_start', 'mag_start', 'length_1', 'identity', 'ref_contig', 'other_contig', 'ref_end','mag_end', 'length_2']
    nucmer_filtered_list = []

    for index, alignment in nucmer_array.iterrows():
        ref_contig = alignment['ref_contig']
        ref_contig = str(ref_contig[1:])
        ref_contig_length = contig_dict[ref_contig]['Length']
        ref_contig_length = int(ref_contig_length)

        other_contig = alignment['other_contig']
        other_contig_length = contig_dict[other_contig]['Length']
        other_contig_length = int(other_contig_length)

        #remove non-end matches
        if int(alignment['mag_start']) >=200 and int(alignment['ref_end']) <= ref_contig_length-200 or int(alignment['ref_start']) >= 200 and int(alignment['mag_end']) <= other_contig_length -200:
            #if int(alignment['mag_start']) < int(alignment['mag_end']) and int(alignment['ref_start']) < int(alignment['ref_end']):
            my_list = [alignment.ref_start, alignment.mag_start, alignment.length_1, alignment.identity, alignment.ref_contig, alignment.other_contig, alignment.ref_end, alignment.mag_end, alignment.length_2]
            nucmer_filtered_list.append(my_list)

    new_array = pd.DataFrame(nucmer_filtered_list, columns=col_names)
    return new_array

