from Bio import SeqIO

def consensus_maker(final_nucmer, mag, rep_mag, x):

    #make nucmer outputs into dictionary with just contig names and what they align to
    global c1, c2, c1s1, c1e1, c1s2, c1e2, c2s, c2e, c3s, c3e, c1_1
    ref_list = final_nucmer['ref_contig'].to_list()
    other_list = final_nucmer['other_contig'].to_list()
    #makes dictionary that looks like {'other': ['ref_1', 'ref_2']} == {'NODE_2767_length_14969_cov_5.423897': ['NODE_1302_length_24101_cov_6.015013', 'NODE_1392_length_23076_cov_6.991834']}
    dict = {}

    for key in other_list:
        for value in ref_list:
            if key not in dict:
                dict[key]=[value[1:]]
                ref_list.remove(value)
                break
            else:
                dict[key].append(value[1:])
                break

    #this list has all of the 'other' contig names but without the duplicates
    #will use for referring to dictionary
    new_other_list = []
    for key in other_list:
        if key not in new_other_list:
                new_other_list.append(key)
    print('Number of keys in dictionary: ',len(new_other_list))
    print('This is the contig match dictionary:', dict)
    print('This is the new list of other contigs:', new_other_list)

    #mag has 'other_contig's and rep_mag has 'ref_contig's
    #names of fasta files
    ref_input = rep_mag + '.' + x
    other_input = mag + '.' + x

    #going through dictionary one key at a time
    for key in new_other_list:
        for index, alignment in final_nucmer.iterrows():
        #get the starts and ends for splicing the sequences
            if key == alignment['other_contig'] and dict.get(key)[0] == alignment['ref_contig'][1:]: #this isnt finding anything
                c1s1 = int(alignment['mag_start']) #contig 1 start 1 (other contig)
                c1e1 = int(alignment['mag_end']) #contig 1 end 1
                c2s = int(alignment['ref_start']) #contig 2 start (first ref contig)
                c2e = int(alignment['ref_end']) #contig 2 end
                for other_fasta in SeqIO.parse(other_input, 'fasta'):
                    if other_fasta.id == key:
                        c1 = other_fasta.seq  # sequence of the key 'other_contig'
                        c1_1 = c1[c1s1:c1e1]  # section of contig aligning with ref 1
                for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
                    if rep_fasta.id == dict.get(key)[0]:
                        c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                        c2 = c2[c2s:c2e]  # section aligning with 'other' contig
            if key == alignment['other_contig'] and dict.get(key)[1] == alignment['ref_contig'][1:]:
                c1s2 = alignment['mag_start'] #contig 1 start 2 (other contig)
                c1e2 = alignment['mag_end'] #contig 1 end 2
                c3s = alignment['ref_start'] #contig 3 end
                c3e = alignment['ref_end'] #contig 3 end
                for other_fasta in SeqIO.parse(other_input, 'fasta'):
                    if other_fasta.id == key:
                        c1 = other_fasta.seq  # sequence of the key 'other_contig'
                        c1_2 = c1[c1s2:c1e2]  # section of contig aligning with ref 2
                for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
                    if rep_fasta.id == dict.get(key)[1]:
                        c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                        c3 = c3[c3s:c3e]  # section aligning with 'other' contig

        #alignment for each key thing and values set in the dictionary
        from Bio import Align
        from Bio.Align import AlignInfo
        from Bio import AlignIO
        aligner = Align.PairwiseAligner()
        alignment_1 = aligner.align(c1_1, c2)
        print(key, ' and ',dict.get(key)[0], ' alignment:', alignment_1)
        alignment_1 = AlignIO.read(alignment_1)
        align1_info = AlignInfo.SummaryInfo(alignment_1)
        consensus_1 = align1_info.dumb_consensus()
        print(key, ' and ',dict.get(key)[0], ' consensus:', consensus_1)

    #
    #         #works to return sequence of named contig
    #         for other_fasta in SeqIO.parse(other_input, 'fasta'):
    #             if other_fasta.id == key:
    #                 c1 = other_fasta.seq #sequence of the key 'other_contig'
    #                 c1_1 = c1[c1s1:c1e1] #section of contig aligning with ref 1
    #                 c1_2 = c1[c1s2:c1e2] #section of contig aligning with ref 2
    #                 print(c1, '\n', c1_1, '\n', c1_2)
    #         for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
    #             if rep_fasta.id == dict.get(key)[0]:
    #                 c2 = rep_fasta.seq #sequence of first 'rep_contig' in value list
    #                 c2 = c2[c2s:c2e] #section aligning with 'other' contig
    #             if rep_fasta.id == dict.get(key)[1]:
    #                 c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
    #                 c3 = c3[c3s:c3e] #section aligning with 'other' contig
    # #


    #         #align
    #         from Bio import Align
    #         aligner = Align.PairwiseAligner()
    #         alignment1 = aligner.align(c1_1, c2)
    #         print('Alignment of other and ref1:', alignment1)

            #
            #      for index, alignment in final_nucmer.iterrows():
            # if key == alignment['other_contig']:
            #     print('dict key', key, ' and other contig ', alignment['other_contig'], 'match' )
            # if dict.get(key)[0] == alignment['ref_contig'][1:]:
            #     print('first ref alignment from dictionary ', dict.get(key)[0], ' and ', alignment['ref_contig'], 'match')

    # if key == alignment['other_contig'] and dict.get(key)[0] == alignment['ref_contig'][
    #                                                             1:]:  # this isnt finding anything
    #     c1s1'] = int(alignment['mag_start'])  # contig 1 start 1 (other contig)
    #     c1e1'] = int(alignment['mag_end'])  # contig 1 end 1
    #     c2s'] = int(alignment['ref_start'])  # contig 2 start (first ref contig)
    #     c2e'] = int(alignment['ref_end'])  # contig 2 end
    # if key == alignment['other_contig'] and dict.get(key)[1] == alignment['ref_contig'][1:]:
    #     c1s2'] = alignment['mag_start']  # contig 1 start 2 (other contig)
    #     c1e2'] = alignment['mag_end']  # contig 1 end 2
    #     c3s'] = alignment['ref_start']  # contig 3 end
    #     c3e'] = alignment['ref_end']  # contig 3 end

    # from Bio import pairwise2
    # from Bio.Align import AlignInfo
    # alignment_1 = pairwise2.align.globalms(c1_1, c2, 2, -1, -.5, -.1, one_alignment_only=True)
    # print(key, ' and ', dict.get(key)[0], ' alignment:', alignment_1)
    # align1_info = AlignInfo.SummaryInfo(alignment_1)
    # consensus_1 = align1_info.dumb_consensus()
    # print(key, ' and ', dict.get(key)[0], ' consensus:', consensus_1)