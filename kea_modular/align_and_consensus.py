from Bio import SeqIO

def consensus_maker(final_nucmer, mag, rep_mag, x):

    #make nucmer outputs into dictionary with just contig names and what they align to
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

    new_other_list = []
    for key in other_list:
        if key not in new_other_list:
                new_other_list.append(key)

    print('This is the contig match dictionary:', dict)
    print('This is the new list of other contigs:', new_other_list)

    #mag has 'other_contig's and rep_mag has 'ref_contig's
    #names of fasta files
    ref_input = rep_mag + '.' + x
    other_input = mag + '.' + x

    for key in new_other_list:
        #works to return sequence of named contig
        for other_fasta in SeqIO.parse(other_input, 'fasta'):
            if other_fasta.id == key:
                c1 = other_fasta.seq #sequence of the key 'other_contig'
                print('Other contig ', key, ' sequence:', repr(c1))
        for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
            if rep_fasta.id == dict.get(str(key))[0]:
                c2 = rep_fasta.seq #sequence of first 'rep_contig' in value list
                print('First ref contig match ', dict.get(key)[0], ' sequence:', repr(c2))
            if rep_fasta.id == dict.get(key)[1]:
                c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
                print('Second ref contig', dict.get(key)[1], ' match sequence:', repr(c3))


# practice = {'NODE_2767_length_14969_cov_5.423897': ['NODE_1302_length_24101_cov_6.015013', 'NODE_1392_length_23076_cov_6.991834']}
# print(practice.get('NODE_2767_length_14969_cov_5.423897')[0])

# for key, value in dict:
#     #works to return sequence of named contig
#     for other_fasta in SeqIO.parse(other_input, 'fasta'):
#         if other_fasta.id == key:
#             c1 = other_fasta.seq #sequence of the key 'other_contig'
#             print('Other contig sequence:', c1)
#     for rep_fasta in SeqIO.parse(ref_input, 'fasta'):
#         if rep_fasta.id == value[0]:
#             c2 = rep_fasta.seq #sequence of first 'rep_contig' in value list
#             print('First ref contig match sequence:', c2)
#         if rep_fasta.id == value[1]:
#             c3 = rep_fasta.seq  # sequence of first 'rep_contig' in value list
#             print('Second ref contig match sequence:', c3)



# #works to return sequence of named contig
# for seq_record in SeqIO.parse('2014_S35.10.fa', 'fasta'):
#     if seq_record.id == 'NODE_22088_length_2502_cov_4.060074':
#         print(repr(seq_record.seq))