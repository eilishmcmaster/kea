from Bio import SeqIO
import os

def consensus_maker(final_nucmer, mag, rep_mag, x):

    #make nucmer outputs into dictionary with just contig names and what they align to
    global c1, c2, c1s1, c1e1, c1s2, c1e2, c2s, c2e, c3s, c3e, c1_1
    bad_ref_list = final_nucmer['ref_contig'].to_list()
    ref_list = []
    for item in bad_ref_list:
        ref_list.append(item[1:])
    death_list = list(ref_list)

    other_list = final_nucmer['other_contig'].to_list()
    #makes dictionary that looks like {'other': ['ref_1', 'ref_2']} == {'NODE_2767_length_14969_cov_5.423897': ['NODE_1302_length_24101_cov_6.015013', 'NODE_1392_length_23076_cov_6.991834']}
    dict = {}
    new_contig_names = []
    new_contig_sequences = []

    for key in other_list:
        for value in death_list:
            if key not in dict:
                dict[key]=[value]
                death_list.remove(value)
                break
            else:
                dict[key].append(value)
                break

    print('dict: ', dict)
    #this list has all of the 'other' contig names but without the duplicates
    #will use for referring to dictionary
    new_other_list = []
    for key in other_list:
        if key not in new_other_list:
                new_other_list.append(key)

    #mag has 'other_contig's and rep_mag has 'ref_contig's
    #names of fasta files
    ref_input = rep_mag + '.' + x
    other_input = mag + '.' + x


#       contig 2 (ref1) (value[0])                                         contig 3 (ref2) (value[1])
# ----------------------------|-------|                            |________|____________________________
#                             |       |                            |        |
#                             |_______|____________________________|________|
#                                       contig 1 (other) (key)


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
                        u2 = c2[:c2s] #unused part of contig 2
            if key == alignment['other_contig'] and dict.get(key)[1] == alignment['ref_contig'][1:]:
                c1s2 = alignment['mag_start'] #contig 1 start 2 (other contig)
                c1e2 = alignment['mag_end'] #contig 1 end 2
                c3s = alignment['ref_start'] #contig 3 end
                c3e = alignment['ref_end'] #contig 3 end
                for other_fasta in SeqIO.parse(other_input, 'fasta'):
                    if other_fasta.id == key:
                        c1 = other_fasta.seq  # sequence of the key 'other_contig'
                        c1_2 = c1[c1s2:c1e2]  # section of contig aligning with ref 2
                        u1 = c1[c1e1:c1s2] #unused section of contig 1
                for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
                    if rep_fasta.id == dict.get(key)[1]:
                        c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                        c3 = c3[c3s:c3e]  # section aligning with 'other' contig
                        u3 = c3[c3e:] #unused section of contig 3

        #alignment for each key thing and values set in the dictionary
        consensus_1 = ''
        for i in range(len(c2)):
            if c1_1[i] == c2[i]:
                consensus_1 = consensus_1 + c1_1[i]
            else:
                consensus_1 = consensus_1 + 'N'

        consensus_2 = ''
        for i in range(len(c3)):
            if c1_2[i] == c3[i]:
                consensus_2 = consensus_2 + c1_2[i]
            else:
                consensus_2 = consensus_2 + 'N'

        #making the final sequence and adding it to the list
        final = u2 + consensus_1 + u1 + consensus_2 + u3
        new_contig_sequences.append(final)
        #making the final sequence name and adding it to the list
        new_contig_name = str(rep_mag +'_gapclosed_' + str(len(final)))
        new_contig_names.append(new_contig_name)

    #EXIT THE DICTIONARY LOOP

    # make a list of the contigs to keep (contigs from the original ref mag without the used contig
    all_contigs_list = []
    with open(ref_input, 'r') as original_fasta:
        for line in original_fasta:
            if '>' in line:
                all_contigs_list.append(line[1:-1])
    unmodified_contig_names = list(set(all_contigs_list) - set(ref_list))  # this is a list with the contigs to write to the new fasta
    unmodified_contig_sequences = []

    #writing sequences to a list
    import re
    string = open(ref_input, 'r').read()
    for i in range(len(unmodified_contig_names)):
        m = re.findall('>' + unmodified_contig_names[i] + '\n(.*?)>', string, re.DOTALL)
        m = str(m).strip("[]").strip('\'')
        unmodified_contig_sequences.append(m)

    #make new fasta and write to it
    new_fasta_name = 'improved_'+ rep_mag + '_using_' + mag #makes improved_rep_using_other
    with open(new_fasta_name, 'w') as new_fasta:
        for i in range(len(new_contig_names)):
            new_fasta.write('>' + new_contig_names[i] + '\n' + new_contig_sequences[i] + '\n')
        for i in range(len(unmodified_contig_names)):
            new_fasta.write('>' + unmodified_contig_names[i]+ '\n' + unmodified_contig_sequences[i])  # write unused contigs -- all contigs in the original reference minus those in the ref_list

    return new_fasta_name #this is a new fasta for each mag x ref improvement -- must be assigned to ref_mag position


