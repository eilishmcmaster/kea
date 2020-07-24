

def gap_overlap(nucmer_array,contig_dict):
    #nucmer array has all the nucmer alignments for the rep mag and other mag in the cluster with low length and identity removed
    #contig dict has all of the contig lengths of all input files in it
    for alignment in nucmer_array:
        if nucmer_array.mag_start <= 200:
            print(alignment)