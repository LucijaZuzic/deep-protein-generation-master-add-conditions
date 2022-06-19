# -*- coding: utf-8 -*-
"""
Created on Sun May 29 10:04:06 2022

@author: Lucija
"""
import matplotlib.pyplot as plt
import numpy as np

luxA_reference = "MKFGNFLLTYQPPQFSQTEVMKRLVKLGRISEECGFDTVWLLEHHFTEFGLLGNPYVAAAYLLGATKKLNVGTAAIVLPTAHPVRQLEDVNLLDQMSKGRFRFGICRGLYNKDFRVFGTDMNNSRALAECWYGLIKNGMTEGYMEADNEHIKFHKVKVNPAAYSRGGAPVYVVAESASTTEWAAQFGLPMILSWIINTNEKKAQLELYNEVAQEYGHDIHNIDHCLSYITSVDHDSIKAKEICRKFLGHWYDSYVNATTIFDDSDQTRGYDFNKGQWRDFVLKGHKDTNRRIDYSYEINPVGTPQECIDIIQKDIDATGISNICCGFEANGTVDEIIASMKLFQSDVMPFLKEKQRSLLY"

def sequence_identity(filename,title):
    try:
        file = open("..\\sequences\\" + filename, 'r')
    except:
        print(filename + "not found")
        return
    lines = file.readlines() 
    file.close()
    print(title)
    lines = [line.strip() for line in lines] 
    difference = []
    score = []
    luxA = luxA_reference + ""
    for sequence_index in range(1, len(lines), 2):
        if (lines[sequence_index - 1] == '>P19839') or (lines[sequence_index - 1] == '>luxAReference') or (lines[sequence_index].replace('-','') == luxA_reference):
            print("found")
            luxA = lines[sequence_index]
            break 
    print(luxA)
    for sequence_index in range(1, len(lines), 2):
        if (lines[sequence_index - 1] == '>P19839') or (lines[sequence_index - 1] == '>luxAReference') or (lines[sequence_index].replace('-','') == luxA_reference):
            print("found") 
            continue 
        difference_count = sum(a != b for a, b in zip(lines[sequence_index], luxA)) + abs(len(luxA) - len(lines[sequence_index]))
        identity = 1 - difference_count / max(len(luxA), len(lines[sequence_index]))
        difference.append(difference_count)
        score.append(identity)
    plt.hist(score)
    plt.xlabel("Percentage of sequence match")
    plt.ylabel("Number of sequences")  
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\sequence_identity\\" + new_filename + "_sequence_identity.png", bbox_inches='tight')
    print(np.min(score), np.max(score), np.mean(score), np.std(score)) 
    plt.close()
    plt.hist(difference)
    plt.xlabel("Number of mutations")
    plt.ylabel("Number of sequences")  
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\difference_count\\" + new_filename + "_difference_count.png", bbox_inches='tight')
    plt.close()
    print(np.min(difference), np.max(difference), np.mean(difference), np.std(difference)) 
    retval = ["",""]
    retval[0] = title + ";" + str(np.round(np.min(score) * 100,3)) + ";" + str(np.round(np.max(score) * 100,3)) + ";" + str(np.round(np.mean(score) * 100,3)) + ";" + str(np.round(np.std(score) * 100,3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(difference),3)) + ";" + str(np.round(np.max(difference),3)) + ";" + str(np.round(np.mean(difference),3)) + ";" + str(np.round(np.std(difference),3))  + "\n" 
    return retval 
 
    
def sequence_identity_multiple(filenames,title):
    lines = []
    all_names = ""
    for filename in filenames:
        try:
            file = open("..\\sequences\\" + filename, 'r')
        except:
            print(filename + "not found")
            return
        lines_part = file.readlines()
        all_names += filename.replace('_ORIGINAL.txt', '_').replace('_ORIGINAL.txt', '_').replace('.txt', '_').replace('.fa', '_').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
        lines_part = [line.strip().replace('-', '') for line in lines_part]
        lines += lines_part 
        file.close()
    print(title) 
    lines = [line.strip() for line in lines] 
    difference = []
    score = []
    luxA = luxA_reference + ""
    for sequence_index in range(1, len(lines), 2):
        if (lines[sequence_index - 1] == '>P19839') or (lines[sequence_index - 1] == '>luxAReference') or (lines[sequence_index].replace('-','') == luxA_reference):
            print("found")
            luxA = lines[sequence_index]
            break 
    print(luxA)
    for sequence_index in range(1, len(lines), 2):
        if (lines[sequence_index - 1] == '>P19839') or (lines[sequence_index - 1] == '>luxAReference') or (lines[sequence_index].replace('-','') == luxA_reference):
            print("found") 
            continue 
        difference_count = sum(a != b for a, b in zip(lines[sequence_index], luxA)) + abs(len(luxA) - len(lines[sequence_index]))
        identity = 1 - difference_count / max(len(luxA), len(lines[sequence_index]))
        difference.append(difference_count)
        score.append(identity)
    plt.hist(score)
    plt.xlabel("Percentage of sequence match")
    plt.ylabel("Number of sequences")  
    plt.title(title)
    plt.savefig("..\\results\\sequence_identity\\" + all_names + "sequence_identity.png", bbox_inches='tight')
    print(np.min(score), np.max(score), np.mean(score), np.std(score)) 
    plt.close()
    plt.hist(difference)
    plt.xlabel("Number of mutations")
    plt.ylabel("Number of sequences")  
    plt.title(title)
    plt.savefig("..\\results\\difference_count\\" + all_names + "difference_count.png", bbox_inches='tight')
    plt.close()
    print(np.min(difference), np.max(difference), np.mean(difference), np.std(difference))
    retval = ["",""]
    retval[0] = title + ";" + str(np.round(np.min(score) * 100,3)) + ";" + str(np.round(np.max(score) * 100,3)) + ";" + str(np.round(np.mean(score) * 100,3)) + ";" + str(np.round(np.std(score) * 100,3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(difference),3)) + ";" + str(np.round(np.max(difference),3)) + ";" + str(np.round(np.mean(difference),3)) + ";" + str(np.round(np.std(difference),3))  + "\n" 
    return retval 
  
output_string_score = "Sequence set;Maximal sequence match (percentage);Maximal sequence match (percentage);Average sequence match (percentage);Standard deviation\n"
output_string_difference = "Sequence set;Minimal number of mutations;Maximal number of mutations;Average number of mutations;Standard deviation\n"
 
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_variants_ORIGINAL.txt','Variants for basic MSA-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol0_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol1_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol2_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_samples_ORIGINAL.txt','Samples for basic MSA-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol0_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol1_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_sol2_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_variants_ORIGINAL.txt','Variants for basic AR-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol0_variants_ORIGINAL.txt','Variants for AR-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol1_variants_ORIGINAL.txt','Variants for AR-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol2_variants_ORIGINAL.txt','Variants for AR-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_samples_ORIGINAL.txt','Samples for basic AR-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol0_samples_ORIGINAL.txt','Samples for AR-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol1_samples_ORIGINAL.txt','Samples for AR-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_sol2_samples_ORIGINAL.txt','Samples for AR-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_with_conditions_sol0_variants_ORIGINAL.txt','Variants for conditional MSA-VAE model with low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_with_conditions_sol1_variants_ORIGINAL.txt','Variants for conditional MSA-VAE model with mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_msavae_with_conditions_sol2_variants_ORIGINAL.txt','Variants for conditional MSA-VAE model with high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_with_conditions_sol0_variants_ORIGINAL.txt','Variants for conditional AR-VAE model with low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_with_conditions_sol1_variants_ORIGINAL.txt','Variants for conditional AR-VAE model with mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('lines_merged\\lines_merged_arvae_with_conditions_sol2_variants_ORIGINAL.txt','Variants for conditional AR-VAE model with high solubility')
output_string_score += s1
output_string_difference += s2
   
s1, s2 = sequence_identity('lines_merged\\lines_merged_PF00296_full.txt', 'Profile of luciferase family')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('training_validation\\luxafilt_llmsa_val.fa', 'Validation data')
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity('training_validation\\luxafilt_llmsa_train.fa', 'Training data')  
output_string_score += s1
output_string_difference += s2
s1, s2 = sequence_identity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], 'Training and validation data')
output_string_score += s1
output_string_difference += s2

output_string_score = output_string_score.replace(".", ",")
file_csv = open("..\\results\\tables\\sequence_identity.csv", "w")
file_csv.write(output_string_score)
file_csv.close()

output_string_difference = output_string_difference.replace(".", ",")
file_csv = open("..\\results\\tables\\difference_count.csv", "w")
file_csv.write(output_string_difference)
file_csv.close()