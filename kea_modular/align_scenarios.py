from Bio import SeqIO

#.......................................................................................................................

#       contig 2 (ref1) (value[0])                                         contig 3 (ref2) (value[1])
# ----------------------------|-------|                            |________|____________________________
#                             |       |                            |        |
#                             |_______|____________________________|________|
#                                       contig 1 (other) (key)

def scenario_1(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict):
    for other_fasta in SeqIO.parse(other_input, 'fasta'):
        if other_fasta.id == key:
            c1 = other_fasta.seq  # sequence of the key 'other_contig'
            c1_1 = c1[c1s1:c1e1]  # section of contig aligning with ref 1
            c1_2 = c1[c1s2:c1e2]  # section of contig aligning with ref 2
            u1 = c1[c1e1:c1s2]  # unused section of contig 1

    for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
        if rep_fasta.id == dict.get(key)[1]:
            c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c3_1 = c3[c3s:c3e]  # section aligning with 'other' contig
            u3 = c3[c3e:]  # unused section of contig 3
        if rep_fasta.id == dict.get(key)[0]:
            c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c2_1 = c2[c2s:c2e]  # section aligning with 'other' contig
            u2 = c2[:c2s]  # unused part of contig 2

    # alignment for each key thing and values set in the dictionary
    consensus_1 = ''
    for i in range(len(c2_1)):
        if c1_1[i] == c2_1[i]:
            consensus_1 = consensus_1 + c1_1[i]
        else:
            consensus_1 = consensus_1 + 'N'

    consensus_2 = ''
    for i in range(len(c3_1)):
        if c1_2[i] == c3_1[i]:
            consensus_2 = consensus_2 + c1_2[i]
        else:
            consensus_2 = consensus_2 + 'N'

    # making the final sequence and adding it to the list
    final = u2 + consensus_1 + u1 + consensus_2 + u3

    return final

#.......................................................................................................................

#         contig 2 (ref1) (value[0])
# -------------------------|---------|--------|      contig 3 (ref2) (value[1])
#                          |         |________|_________________________________
#                          |_________|________|____________|
#                               contig 1 (other) (key)

def scenario_2(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict):
    for other_fasta in SeqIO.parse(other_input, 'fasta'):
        if other_fasta.id == key:
            c1 = other_fasta.seq  # sequence of the key 'other_contig'
            c1_1 = c1[c1s1:c1s2]  # section of contig aligning with ref 1
            c1_2 = c1[c1s2:c1e1]  # section of contig aligning with ref 2
            c1_3 = c1[c1e1:c1e2]  # unused section of contig 1

    for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
        if rep_fasta.id == dict.get(key)[1]:
            c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c3_1 = c3[c3s:(c3s+c1e1-c1s2)]  # section aligning with 'other' contig
            c3_2 = c3[(c3s+c1e1-c1s2):c3e]  # unused section of contig 3
            u2 = c3[c3s:]
        if rep_fasta.id == dict.get(key)[0]:
            c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c2_1 = c2[c2s:(c2s+c1s2-c1s1)]  # section aligning with 'other' contig
            c2_2 = c2[(c2s+c1s2-c1s1):c2e]  # unused part of contig 2
            u1 = c2[:c2s]
    print('c1_1, c2_1: ',len(c1_1), len(c2_1))
    print('middle: ', len(c2_2), len(c3_1), len(c1_2))
    print('c22: ', c2_2[:10], c2_2[-10:])
    print('c31: ', c3_1[:10], c3_1[-10:])
    print('c12: ', c1_2[:10], c1_2[-10:])
    print('end; ', len(c3_2), len(c1_3))
    print('how long total should be: ', c1e2-c1s1)

    consensus_1 = ''
    for i in range(len(c2_1)):
        if c1_1[i] == c2_1[i]:
            consensus_1 = consensus_1 + c1_1[i]
        else:
            consensus_1 = consensus_1 + 'N'
    consensus_3 = ''
    for i in range(len(c3_2)):
        if c1_3[i] == c3_2[i]:
            consensus_3 = consensus_3 + c1_3[i]
        else:
            consensus_3 = consensus_3 + 'N'

    consensus_2 = ''
    for i in range(len(c2_2)):
        if c2_2[i] == c3_1[i] == c1_2[i]: #all things match
            consensus_2 = consensus_2 + c1_2[i]
        elif c2_2[i] == c3_1[i] != c1_2[i]: #two match
            consensus_2 = consensus_2 + c2_2[i]
        elif c2_2[i] != c3_1[i] == c1_2[i]: #two match
            consensus_2 = consensus_2 + c1_2[i]
        elif c2_2[i] == c1_2[i] != c3_1[i] :  #two match
            consensus_2 = consensus_2 + c2_2[i]
        else:
            consensus_2 = consensus_2 + 'N'

    final = u1 + consensus_1 + consensus_2 + consensus_3 + u2
    return final