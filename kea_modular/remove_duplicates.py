import numpy as np
def remove_dupes(filtered_nucmer):
    # make a pivot table with the number of duplicates in the ref_contig column
    dupes = filtered_nucmer.pivot_table(index=['ref_contig'], aggfunc='size')

    # retutrn to a dataframe and set the count column 0 as integers
    df = dupes.reset_index()
    df[0] = df[0].astype(int)
    # make a list of columns to remove that have more than 2 duplicate values
    contigs_to_remove = []
    for index, info in df.iterrows():
        if info[0] >= 2:
            ref_contig = info['ref_contig']
            ref_contig = str(' '+ ref_contig[1:])
            contigs_to_remove.append(ref_contig)
    # remove the >2 duplicates form the filtered_nucmer
    filtered_nucmer = filtered_nucmer[~filtered_nucmer['ref_contig'].isin(contigs_to_remove)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # make a pivot table with the number of duplicates in the ref_contig column
    dupes2 = filtered_nucmer.pivot_table(index=['other_contig'], aggfunc='size')

    # retutrn to a dataframe and set the count column 0 as integers
    df2 = dupes2.reset_index()
    df2[0] = df2[0].astype(int)
    contigs_to_remove2 = []
    for index, info in df2.iterrows():
        if info[0] != 2:
            other_contig2 = info['other_contig']
            other_contig2 = str(other_contig2)
            contigs_to_remove2.append(other_contig2)

    # remove the >2 duplicates form the filtered_nucmer
    filtered_nucmer = filtered_nucmer[~filtered_nucmer['other_contig'].isin(contigs_to_remove2)]
    #sort so that the dictionary in the next function will work right
    filtered_nucmer = filtered_nucmer.sort_values('ref_start', ascending=False)

    print(filtered_nucmer)
    start_end_dict = {}
    for index, alignment in filtered_nucmer.iterrows():
        if int(alignment['mag_start'])>=200:
            start_end_dict.update({str(alignment['ref_contig']):[0,str(alignment['other_contig'])]}) #match is start of ref contig =0
        else:
            start_end_dict.update({str(alignment['ref_contig']):[1,str(alignment['other_contig'])]}) #match is at end of ref contig = 1

    print('start end dict: ', start_end_dict)

    # make dictionary that has {'[1, other contig]':[ref1, ref2]}
    flipped = {}
    for key, value in start_end_dict.items():
        if str(value) not in flipped:
            flipped[str(value)] = [key]
        else:
            flipped[str(value)].append(key)

    length_dict = {key: len(value) for key, value in flipped.items()}
    other_remove_list = []  # list of other contig things to remove
    ref_remove_list = []  # list of reference contigs to remove
    # find other contig amtches and extract
    for key, value in length_dict.items():
        if value >= 2:
            other_remove_list.append(key[5:-2])
    # make a list of the reference contigs to remove from array
    for item in other_remove_list:
        for key, value in start_end_dict.items():
            if item == value[1]:
                ref_remove_list.append(key)

    filtered_nucmer = filtered_nucmer[~filtered_nucmer['ref_contig'].isin(ref_remove_list)]
    print(filtered_nucmer)

    return filtered_nucmer
