# -*- coding: utf-8 -*-
"""
Created on Sun May 22 16:50:29 2022

@author: Lucija
"""
import matplotlib.pyplot as plt
import numpy as np

luxA_reference = "MKFGNFLLTYQPPQFSQTEVMKRLVKLGRISEECGFDTVWLLEHHFTEFGLLGNPYVAAAYLLGATKKLNVGTAAIVLPTAHPVRQLEDVNLLDQMSKGRFRFGICRGLYNKDFRVFGTDMNNSRALAECWYGLIKNGMTEGYMEADNEHIKFHKVKVNPAAYSRGGAPVYVVAESASTTEWAAQFGLPMILSWIINTNEKKAQLELYNEVAQEYGHDIHNIDHCLSYITSVDHDSIKAKEICRKFLGHWYDSYVNATTIFDDSDQTRGYDFNKGQWRDFVLKGHKDTNRRIDYSYEINPVGTPQECIDIIQKDIDATGISNICCGFEANGTVDEIIASMKLFQSDVMPFLKEKQRSLLY"
BLOSSUM_MATRIX = []
characters = []
lambda_value = np.log(2) / 2

def read_matrix():
    file = open("BLOSSUM62_matrix.txt")
    lines = file.readlines()  
    for line_num in range(len(lines)):
        new_line = lines[line_num].strip().replace("\n", "")
        while new_line.count("  ") != 0:
            new_line = new_line.replace("  ", " ") 
        new_line = new_line.split(" ")
        if line_num == 0:
            for char in new_line:
                characters.append(char)
        else:
            new_line = new_line[1:]
            new_line = [int(element) for element in new_line] 
            BLOSSUM_MATRIX.append(new_line) 
    characters[-1] = '-'
            
read_matrix()

def similarity(filename,title):
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
        BLOSSUM_score = sum(BLOSSUM_MATRIX[characters.index(a)][characters.index(b)] for a, b in zip(lines[sequence_index], luxA))
        if len(luxA) > len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(luxA[i])][characters.index('-')] for i in range(len(lines[sequence_index]), len(luxA)))
        if len(luxA) < len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(lines[sequence_index][i])][characters.index('-')] for i in range(len(luxA), len(lines[sequence_index])))
        chance = (2 ** (BLOSSUM_score * lambda_value)) / (1 + (2 ** (BLOSSUM_score * lambda_value)))
        difference.append(chance)
        score.append(BLOSSUM_score)
    plt.hist(score)
    plt.xlabel("BLOSSUM 62 score")
    plt.ylabel("Number of sequences")
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\BLOSSUM62\\" + new_filename + "_BLOSSUM62.png", bbox_inches='tight')
    print(np.min(score), np.max(score), np.mean(score), np.std(score)) 
    plt.close()
    plt.hist(difference)
    plt.xlabel("Probability of relation")
    plt.ylabel("Number of sequences")
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\log_odd\\" + new_filename + "_log_odd.png", bbox_inches='tight')
    plt.close()
    print(np.min(difference), np.max(difference), np.mean(difference), np.std(difference)) 
    retval = ["",""]
    retval[0] = title + ";" + str(np.round(np.min(score),3)) + ";" + str(np.round(np.max(score),3)) + ";" + str(np.round(np.mean(score),3)) + ";" + str(np.round(np.std(score),3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(difference) * 100,3)) + ";" + str(np.round(np.max(difference) * 100,3)) + ";" + str(np.round(np.mean(difference) * 100,3)) + ";" + str(np.round(np.std(difference) * 100,3))  + "\n" 
    return retval      
 
    
def similarity_multiple(filenames,title):
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
        BLOSSUM_score = sum(BLOSSUM_MATRIX[characters.index(a)][characters.index(b)] for a, b in zip(lines[sequence_index], luxA))
        if len(luxA) > len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(luxA[i])][characters.index('-')] for i in range(len(lines[sequence_index]), len(luxA)))
        if len(luxA) < len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(lines[sequence_index][i])][characters.index('-')] for i in range(len(luxA), len(lines[sequence_index])))
        chance = (2 ** (BLOSSUM_score * lambda_value)) / (1 + (2 ** (BLOSSUM_score * lambda_value)))
        difference.append(chance)
        score.append(BLOSSUM_score)
    plt.hist(score)
    plt.xlabel("BLOSSUM 62 score")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\BLOSSUM62\\" + all_names + "BLOSSUM62.png", bbox_inches='tight')
    print(np.min(score), np.max(score), np.mean(score), np.std(score)) 
    plt.close()
    plt.hist(difference)
    plt.xlabel("Probability of relation")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\log_odd\\" + all_names + "log_odd.png", bbox_inches='tight')
    plt.close()
    print(np.min(difference), np.max(difference), np.mean(difference), np.std(difference))
    retval = ["",""]
    retval[0] = title + ";" + str(np.round(np.min(score),3)) + ";" + str(np.round(np.max(score),3)) + ";" + str(np.round(np.mean(score),3)) + ";" + str(np.round(np.std(score),3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(difference) * 100,3)) + ";" + str(np.round(np.max(difference) * 100,3)) + ";" + str(np.round(np.mean(difference) * 100,3)) + ";" + str(np.round(np.std(difference) * 100,3))  + "\n" 
    return retval 
  
output_string_score = "Sequence set;Minimal BLOSSUM 62 score;Maximal BLOSSUM 62 rezultat;Average BLOSSUM 62 score;Standard deviation\n"
output_string_difference = "Sequence set;Minimal probability of relation (percentage);Maximal probability of relation (percentage);Average probability of relation (percentage);Standard deviation\n"
 
s1, s2 = similarity('lines_merged\\lines_merged_msavae_variants_ORIGINAL.txt','Variants for basic MSA-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol0_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol1_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol2_variants_ORIGINAL.txt','Variants for MSA-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = similarity('lines_merged\\lines_merged_msavae_samples_ORIGINAL.txt','Samples for basic MSA-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol0_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol1_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_sol2_samples_ORIGINAL.txt','Samples for MSA-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = similarity('lines_merged\\lines_merged_arvae_variants_ORIGINAL.txt','Variants for basic AR-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol0_variants_ORIGINAL.txt','Variants for AR-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol1_variants_ORIGINAL.txt','Variants for AR-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol2_variants_ORIGINAL.txt','Variants for AR-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = similarity('lines_merged\\lines_merged_arvae_samples_ORIGINAL.txt','Samples for basic AR-VAE model')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol0_samples_ORIGINAL.txt','Samples for AR-VAE model trained on low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol1_samples_ORIGINAL.txt','Samples for AR-VAE model trained on mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_sol2_samples_ORIGINAL.txt','Samples for AR-VAE model trained on high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol0_variants_ORIGINAL.txt','Variants for conditional MSA-VAE model with low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol1_variants_ORIGINAL.txt','Variants for conditional MSA-VAE mode with mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol2_variants_ORIGINAL.txt','Variants for conditional MSA-VAE mode with high solubility')
output_string_score += s1
output_string_difference += s2

s1, s2 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol0_variants_ORIGINAL.txt','Variants for conditional AR-VAE model with low solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol1_variants_ORIGINAL.txt','Variants for conditional AR-VAE mode with mid solubility')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol2_variants_ORIGINAL.txt','Variants for conditional AR-VAE mode with high solubility')
output_string_score += s1
output_string_difference += s2
   
s1, s2 = similarity('lines_merged\\lines_merged_PF00296_full.txt', 'Profile of luciferase family') 
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('training_validation\\luxafilt_llmsa_val.fa', 'Validation data')
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity('training_validation\\luxafilt_llmsa_train.fa', 'Training data')  
output_string_score += s1
output_string_difference += s2
s1, s2 = similarity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], 'Training and validation data')
output_string_score += s1
output_string_difference += s2

output_string_score = output_string_score.replace(".", ",")
file_csv = open("..\\results\\tables\\BLOSSUM62.csv", "w")
file_csv.write(output_string_score)
file_csv.close()

output_string_difference = output_string_difference.replace(".", ",")
file_csv = open("..\\results\\tables\\log_odd.csv", "w")
file_csv.write(output_string_difference)
file_csv.close()