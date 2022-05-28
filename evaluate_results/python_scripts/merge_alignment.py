# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def merge_alignment(filename):
    try:
        file = open("..\\sequences\\aligned_sequences\\" + filename, 'r')
    except:
        return
    print(filename)
    lines = file.readlines()
    file.close()
    new_file_contents = ""
    new_seq = ""
    for line in lines:
        if line[0] == '>':
            if new_seq != "":
                new_file_contents += new_seq + "\n"
            new_file_contents += line.replace('\n', '') + "\n"
            new_seq = ""
        else:
            new_seq += line.replace('\n', '')
    if new_seq != "":
        new_file_contents += new_seq + "\n"
    new_seq = ""
    new_file = open("..\\sequences\\lines_merged\\lines_merged_" + filename, 'w')
    new_file.write(new_file_contents)
    new_file.close()

merge_alignment("PF00296_full.txt")
merge_alignment("arvae_variants_ORIGINAL.txt")
merge_alignment("arvae_sol0_variants_ORIGINAL.txt")
merge_alignment("arvae_sol1_variants_ORIGINAL.txt")
merge_alignment("arvae_sol2_variants_ORIGINAL.txt")
merge_alignment("arvae_samples_ORIGINAL.txt")
merge_alignment("arvae_sol0_samples_ORIGINAL.txt")
merge_alignment("arvae_sol1_samples_ORIGINAL.txt")
merge_alignment("arvae_sol2_samples_ORIGINAL.txt") 
merge_alignment("msavae_with_conditions_sol0_variants_ORIGINAL.txt")
merge_alignment("msavae_with_conditions_sol1_variants_ORIGINAL.txt")
merge_alignment("msavae_with_conditions_sol2_variants_ORIGINAL.txt")
merge_alignment("msavae_variants_ORIGINAL.txt")
merge_alignment("msavae_sol0_variants_ORIGINAL.txt")
merge_alignment("msavae_sol1_variants_ORIGINAL.txt")
merge_alignment("msavae_sol2_variants_ORIGINAL.txt")
merge_alignment("msavae_samples_ORIGINAL.txt")
merge_alignment("msavae_sol0_samples_ORIGINAL.txt")
merge_alignment("msavae_sol1_samples_ORIGINAL.txt")
merge_alignment("msavae_sol2_samples_ORIGINAL.txt") 
merge_alignment("msavae_with_conditions_sol0_variants_ORIGINAL.txt")
merge_alignment("msavae_with_conditions_sol1_variants_ORIGINAL.txt")
merge_alignment("msavae_with_conditions_sol2_variants_ORIGINAL.txt")