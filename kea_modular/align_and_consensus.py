
import os
from align_scenarios import scenario_1, scenario_2, datetime

def consensus_maker(final_nucmer, mag, rep_mag, x, contig_dict):

    #make nucmer outputs into dictionary with just contig names and what they align to
    u1, u2, u3, c1, c2, c3, c1s1, c1e1, c1s2, c1e2, c2s, c2e, c3s, c3e, c1_1, c1_2, c2_1, c3_1 = [""]*18
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
            if key == alignment['other_contig'] and dict.get(key)[0] == alignment['ref_contig'][1:]:
                c1s1 = int(alignment['mag_start']) #contig 1 start 1 (other contig)
                c1e1 = int(alignment['mag_end']) #contig 1 end 1
                c2s = int(alignment['ref_start']) #contig 2 start (first ref contig)
                c2e = int(alignment['ref_end']) #contig 2 end
            if key == alignment['other_contig'] and dict.get(key)[1] == alignment['ref_contig'][1:]:
                c1s2 = int(alignment['mag_start']) #contig 1 start 2 (other contig)
                c1e2 = int(alignment['mag_end']) #contig 1 end 2
                c3s = int(alignment['ref_start']) #contig 3 end
                c3e = int(alignment['ref_end']) #contig 3 end

        if c1e1<c1s2: #simple, reference contigs dont overlap
            print(datetime.now(), 'Simple alignment initiated (scenario 1)')
            final = scenario_1(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict)
        elif c1e1>c1s2: #references overlap
            print(datetime.now(), 'Overlapping contig alignment initiated (scenario 2)')
            final = scenario_2(other_input, ref_input, c1s1,c1e1,c1s2,c1e2,c3s,c3e,c2s,c2e,key,dict)

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

    unmodified_contig_names = [x for x in all_contigs_list if x not in ref_list] # this is a list with the contigs to write to the new fasta
    unmodified_contig_sequences = []

    #writing sequences to a list
    import re
    string = open(ref_input, 'r').read()
    for i in range(len(unmodified_contig_names)):
        m = re.findall('>' + unmodified_contig_names[i] + '\n(.*?)\n>', string, re.DOTALL)
        m = str(m).strip("[]").strip('\'')
        unmodified_contig_sequences.append(m)
    rec = string.split(unmodified_contig_names[-1])[1] #get the last sequence in the file
    unmodified_contig_sequences.append(rec)
    unmodified_contig_sequences = [s.replace('\\n', '') for s in unmodified_contig_sequences]
    unmodified_contig_sequences = [s.replace('\n', '') for s in unmodified_contig_sequences]

    #make new fasta and write to it
    new_fasta_name = 'improved_'+ rep_mag  #makes improved_rep_#

    new_fasta_name_x = new_fasta_name + '.' + x
    with open(new_fasta_name_x, 'w') as new_fasta:
        for i in range(len(new_contig_names)):
            new_fasta.write(str('>' + new_contig_names[i] + '\n' + new_contig_sequences[i] + '\n'))
        for i in range(len(unmodified_contig_names)):
            new_fasta.write('>' + unmodified_contig_names[i]+ '\n' + unmodified_contig_sequences[i] + '\n')  # write unused contigs -- all contigs in the original reference minus those in the ref_list


    #adding new contig lengths to contig_dict
    improved_contig_dict = {}
    for cont in new_contig_names:
        index = new_contig_names.index(cont)
        length = len(new_contig_sequences[index])
        improved_contig_dict[cont]={'Length':length}

    from deepmerge import always_merger
    contig_dict = always_merger.merge(contig_dict, improved_contig_dict)


    return new_fasta_name, contig_dict #this is a new fasta for each mag x ref improvement -- must be assigned to ref_mag position


