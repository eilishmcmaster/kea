from Bio import SeqIO
from datetime import datetime
#angerey
#.......................................................................................................................

#       contig 2 (ref1) (value[0])                                         contig 3 (ref2) (value[1])
# ----------------------------|-------|                            |________|____________________________
#                             |       |                            |        |
#                             |_______|____________________________|________|
#                                       contig 1 (other) (key)

def scenario_1(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict):
    contig1 = ''
    contig2 = ''
    contig3 = ''
    for other_fasta in SeqIO.parse(other_input, 'fasta'):
        if other_fasta.id == key:
            print(datetime.now(), 'Closing gap using other contig ', key)
            c1 = other_fasta.seq  # sequence of the key 'other_contig'
            c1_1 = c1[c1s1:c1e1]  # section of contig aligning with ref 1
            c1_2 = c1[c1s2:c1e2]  # section of contig aligning with ref 2
            u1 = c1[c1e1:c1s2]  # unused section of contig 1
            contig3 = key

    for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
        if rep_fasta.id == dict.get(key)[1]:
            print(datetime.now(), 'Reference contig 1: ', dict.get(key)[1])
            c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c3_1 = c3[c3s:c3e]  # section aligning with 'other' contig
            u3 = c3[c3e:]  # unused section of contig 3
            contig2 = dict.get(key)[1]
        if rep_fasta.id == dict.get(key)[0]:
            print(datetime.now(), 'Reference contig 2: ', dict.get(key)[0])
            c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
            c2_1 = c2[c2s:c2e]  # section aligning with 'other' contig
            u2 = c2[:c2s]  # unused part of contig 2
            contig1 = dict.get(key)[0]

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

    with open('../log.tsv', 'a') as stupid:
        stupid.write(ref_input + '\t' + contig1 + '\t' + contig2 + '\t'+ other_input +'\t'+ contig3 + '\tscenario_1\t')

    return final


#.......................................................................................................................

#         contig 2 (ref1) (value[0])
# -------------------------|---------|--------|      contig 3 (ref2) (value[1])
#                          |         |________|_________________________________
#                          |_________|________|____________|
#                               contig 1 (other) (key)

def scenario_2(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict):
    contig1 = ''
    contig2 = ''
    contig3 = ''

    # swap_s1=c1s1
    # swap_s2=c1s2
    # swap_e1=c1e1
    # swap_e2=c1e2
    # c1s1=''
    # c1e1=''
    # c1s2=''
    # c1e2=''
    #
    # c1s1=swap_s2
    # c1s2=swap_s1
    # c1e1=swap_e2
    # c1e2=swap_e1
    if c1s2>c1s1:
        for other_fasta in SeqIO.parse(other_input, 'fasta'):
            if other_fasta.id == key:
                print(datetime.now(), 'Closing gap using other contig ', key)
                c1 = other_fasta.seq  # sequence of the key 'other_contig'
                c1_1 = c1[c1s1:c1s2]  # section of contig aligning with ref 1
                c1_2 = c1[c1s2:c1e1]  # section of contig aligning with ref 2
                c1_3 = c1[c1e1:c1e2]  # unused section of contig 1
                contig3 = key

        for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
            if rep_fasta.id == dict.get(key)[1]:
                print(datetime.now(), 'Reference contig 1: ', dict.get(key)[1])
                c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                c3_1 = c3[c3s:(c3s+c1e1-c1s2)]  # section aligning with 'other' contig
                c3_2 = c3[(c3s+c1e1-c1s2):c3e]  # unused section of contig 3
                u2 = c3[c3s:]
                contig2 = dict.get(key)[1]
            if rep_fasta.id == dict.get(key)[0]:
                print(datetime.now(), 'Reference contig 2: ', dict.get(key)[0])
                c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                c2_1 = c2[c2s:(c2s+c1s2-c1s1)]  # section aligning with 'other' contig
                c2_2 = c2[(c2s+c1s2-c1s1):c2e]  # unused part of contig 2
                u1 = c2[:c2s]
                contig1 = dict.get(key)[0]
    else:
        start2=c1s1
        start1=c1s2
        for other_fasta in SeqIO.parse(other_input, 'fasta'):
            if other_fasta.id == key:
                print(datetime.now(), 'Closing gap using other contig ', key)
                c1 = other_fasta.seq  # sequence of the key 'other_contig'
                c1_1 = c1[start2:start1]  # section of contig aligning with ref 1
                c1_2 = c1[start1:c1e1]  # section of contig aligning with ref 2
                c1_3 = c1[c1e1:c1e2]  # unused section of contig 1
                contig3 = key

        for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
            if rep_fasta.id == dict.get(key)[1]:
                print(datetime.now(), 'Reference contig 1: ', dict.get(key)[1])
                c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                c3_1 = c3[c3s:(c3s+c1e1-start1)]  # section aligning with 'other' contig
                c3_2 = c3[(c3s+c1e1-start1):c3e]  # unused section of contig 3
                u2 = c3[c3s:]
                contig2 = dict.get(key)[1]
            if rep_fasta.id == dict.get(key)[0]:
                print(datetime.now(), 'Reference contig 2: ', dict.get(key)[0])
                c2 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                c2_1 = c2[c2s:(c2s+start1-start2)]  # section aligning with 'other' contig
                c2_2 = c2[(c2s+start1-start2):c2e]  # unused part of contig 2
                u1 = c2[:c2s]
                contig1 = dict.get(key)[0]


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
    with open('../log.tsv', 'a') as stupid:
        stupid.write(ref_input + '\t' + contig1 + '\t' + contig2 + '\t'+ other_input +'\t'+ contig3 + '\tscenario_2\t')
    return final