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
    print('df2:', df2)
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

    return filtered_nucmer
