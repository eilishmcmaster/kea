def remove_dupes(filtered_nucmer):
    # make a pivot table with the number of duplicates in the ref_contig column
    dupes = filtered_nucmer.pivot_table(index=['ref_contig'], aggfunc='size')

    # retutrn to a dataframe and set the count column 0 as integers
    df = dupes.reset_index()
    df[0] = df[0].astype(int)

    # make a list of columns to remove that have more than 2 duplicate values
    contigs_to_remove = []
    for index, info in df.iterrows():
        if info[0] > 2:
            ref_contig = info['ref_contig']
            ref_contig = str(ref_contig[1:])
            contigs_to_remove.append(ref_contig)

    # remove the >2 duplicates form the filtered_nucmer
    filtered_nucmer = filtered_nucmer[~filtered_nucmer['ref_contig'].isin(contigs_to_remove)]

    return filtered_nucmer
