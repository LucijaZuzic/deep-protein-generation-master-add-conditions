# -*- coding: utf-8 -*-
"""
Created on Sun May 22 16:50:29 2022

@author: Lucija
"""
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
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

def similarity(filename, cosine_file, cosine_reference, title, extension = ''): 
    try:
        file = open("..\\sequences\\" + filename, 'r')
    except:
        print(filename + "not found")
        return
    lines = file.readlines() 
    file.close()
    print(title)
    lines = [line.strip() for line in lines] 
    score = []
    cosine = []
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
        index_of_cosine = int((sequence_index - 1) / 2)
        if lines[sequence_index - 1][1] == 's' or lines[sequence_index - 1][1] == 'l':
            index_of_cosine = lines[sequence_index - 1].replace('>s', '').replace('>luxa_var', '')
            index_of_cosine = int(index_of_cosine) - 1
        vector = cosine_file[index_of_cosine]
        cosine_similarity = np.dot(vector,cosine_reference)/(np.linalg.norm(vector)*np.linalg.norm(cosine_reference))
        cosine.append(cosine_similarity)
        BLOSSUM_score = sum(BLOSSUM_MATRIX[characters.index(a)][characters.index(b)] for a, b in zip(lines[sequence_index], luxA))
        if len(luxA) > len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(luxA[i])][characters.index('-')] for i in range(len(lines[sequence_index]), len(luxA)))
        if len(luxA) < len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(lines[sequence_index][i])][characters.index('-')] for i in range(len(luxA), len(lines[sequence_index])))
        chance = (2 ** (BLOSSUM_score * lambda_value)) / (1 + (2 ** (BLOSSUM_score * lambda_value)))
        score.append(chance)
    plt.scatter(cosine, score) 
    plt.ylabel("BLOSSUM 62 score")
    plt.xlabel("Cosine similarity")  
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\cosine\\" + new_filename + "_cosine" + extension + ".png", bbox_inches='tight')
    print(np.min(cosine), np.max(cosine), np.mean(cosine), np.std(cosine)) 
    plt.close() 
    retval = title + ";" + str(np.round(np.min(cosine) * 100,3)) + ";" + str(np.round(np.max(cosine) * 100,3)) + ";" + str(np.round(np.mean(cosine) * 100,3)) + ";" + str(np.round(np.std(cosine) * 100,3))  + "\n" 
    return retval      
 
    
def similarity_multiple(filenames, cosine_file, cosine_reference, title, extension=''): 
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
    score = []
    cosine = []
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
        index_of_cosine = int((sequence_index - 1) / 2)
        if lines[sequence_index - 1][1] == 's' or lines[sequence_index - 1][1] == 'l':
            index_of_cosine = lines[sequence_index - 1].replace('>s', '').replace('>luxa_var', '')
            index_of_cosine = int(index_of_cosine) - 1
        vector = cosine_file[index_of_cosine]
        cosine_similarity = np.dot(vector,cosine_reference)/(np.linalg.norm(vector)*np.linalg.norm(cosine_reference))
        cosine.append(cosine_similarity)
        BLOSSUM_score = sum(BLOSSUM_MATRIX[characters.index(a)][characters.index(b)] for a, b in zip(lines[sequence_index], luxA))
        if len(luxA) > len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(luxA[i])][characters.index('-')] for i in range(len(lines[sequence_index]), len(luxA)))
        if len(luxA) < len(lines[sequence_index]):
            BLOSSUM_score += sum(BLOSSUM_MATRIX[characters.index(lines[sequence_index][i])][characters.index('-')] for i in range(len(luxA), len(lines[sequence_index])))
        chance = (2 ** (BLOSSUM_score * lambda_value)) / (1 + (2 ** (BLOSSUM_score * lambda_value)))
        score.append(chance)
    plt.scatter(cosine, score) 
    plt.ylabel("BLOSSUM 62 score")
    plt.xlabel("Cosine similarity")  
    plt.title(title)
    new_filename = filename.replace('_ORIGINAL.txt', '').replace('_ORIGINAL.txt', '').replace('.txt', '').replace('.fa', '').replace('lines_merged\\lines_merged_', '').replace('training_validation\\', '')
    plt.savefig("..\\results\\cosine\\" + new_filename + "_cosine" + extension + ".png", bbox_inches='tight')
    print(np.min(cosine), np.max(cosine), np.mean(cosine), np.std(cosine)) 
    plt.close() 
    retval = title + ";" + str(np.round(np.min(cosine) * 100,3)) + ";" + str(np.round(np.max(cosine) * 100,3)) + ";" + str(np.round(np.mean(cosine) * 100,3)) + ";" + str(np.round(np.std(cosine) * 100,3))  + "\n" 
    return retval     

pca = PCA(n_components=2)

msa_outputs = np.load('../../output/generated_PCA/msa_outputs.npy')[0] 
 
msa_reference = msa_outputs[73]

msa_outputs_samples = np.load('../../output/generated_PCA/msa_outputs_samples.npy')[0]
msa_outputs_variants = np.load('../../output/generated_PCA/msa_outputs_variants.npy')[0]

msa_outputs_sol0 = np.load('../../output/generated_PCA/msa_outputs_sol0.npy')[0]
msa_outputs_sol0_samples = np.load('../../output/generated_PCA/msa_outputs_sol0_samples.npy')[0]
msa_outputs_sol0_variants = np.load('../../output/generated_PCA/msa_outputs_sol0_variants.npy')[0]

msa_outputs_sol1 = np.load('../../output/generated_PCA/msa_outputs_sol1.npy')[0]
msa_outputs_sol1_samples = np.load('../../output/generated_PCA/msa_outputs_sol1_samples.npy')[0]
msa_outputs_sol1_variants = np.load('../../output/generated_PCA/msa_outputs_sol1_variants.npy')[0]

msa_outputs_sol2 = np.load('../../output/generated_PCA/msa_outputs_sol2.npy')[0]
msa_outputs_sol2_samples = np.load('../../output/generated_PCA/msa_outputs_sol2_samples.npy')[0]
msa_outputs_sol2_variants = np.load('../../output/generated_PCA/msa_outputs_sol2_variants.npy')[0]

msa_outputs_cond = np.load('../../output/generated_PCA/msa_outputs_cond.npy')[0]

msa_reference_cond = msa_outputs_cond[73]

msa_outputs_cond_sol0 = np.load('../../output/generated_PCA/msa_outputs_cond_sol0.npy')[0]
msa_outputs_cond_sol1 = np.load('../../output/generated_PCA/msa_outputs_cond_sol1.npy')[0]
msa_outputs_cond_sol2 = np.load('../../output/generated_PCA/msa_outputs_cond_sol2.npy')[0]

ar_outputs = np.load('../../output/generated_PCA/ar_outputs.npy')[0] 

ar_reference = ar_outputs[73]

ar_outputs_samples = np.load('../../output/generated_PCA/ar_outputs_samples.npy')[0]
ar_outputs_variants = np.load('../../output/generated_PCA/ar_outputs_variants.npy')[0]

ar_outputs_sol0 = np.load('../../output/generated_PCA/ar_outputs_sol0.npy')[0]
ar_outputs_sol0_samples = np.load('../../output/generated_PCA/ar_outputs_sol0_samples.npy')[0]
ar_outputs_sol0_variants = np.load('../../output/generated_PCA/ar_outputs_sol0_variants.npy')[0]

ar_outputs_sol1 = np.load('../../output/generated_PCA/ar_outputs_sol1.npy')[0]
ar_outputs_sol1_samples = np.load('../../output/generated_PCA/ar_outputs_sol1_samples.npy')[0]
ar_outputs_sol1_variants = np.load('../../output/generated_PCA/ar_outputs_sol1_variants.npy')[0]

ar_outputs_sol2 = np.load('../../output/generated_PCA/ar_outputs_sol2.npy')[0]
ar_outputs_sol2_samples = np.load('../../output/generated_PCA/ar_outputs_sol2_samples.npy')[0]
ar_outputs_sol2_variants = np.load('../../output/generated_PCA/ar_outputs_sol2_variants.npy')[0]

ar_outputs_cond = np.load('../../output/generated_PCA/ar_outputs_cond.npy')[0]

ar_reference_cond = ar_outputs_cond[73] 

ar_outputs_cond_sol0 = np.load('../../output/generated_PCA/ar_outputs_cond_sol0.npy')[0]
ar_outputs_cond_sol1 = np.load('../../output/generated_PCA/ar_outputs_cond_sol1.npy')[0]
ar_outputs_cond_sol2 = np.load('../../output/generated_PCA/ar_outputs_cond_sol2.npy')[0]

output_string_score = "Sequence set;Minimal cosine similarity;Maximal cosine similarityt;Average cosine similarity;Standard deviation\n"
 
s1 = similarity('lines_merged\\lines_merged_msavae_variants_ORIGINAL.txt', msa_outputs_variants, msa_reference, 'Variants for basic MSA-VAE model')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol0_variants_ORIGINAL.txt', msa_outputs_sol0_variants, msa_reference, 'Variants for MSA-VAE model trained on low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol1_variants_ORIGINAL.txt', msa_outputs_sol1_variants, msa_reference, 'Variants for MSA-VAE model trained on mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol2_variants_ORIGINAL.txt', msa_outputs_sol2_variants, msa_reference, 'Variants for MSA-VAE model trained on high solubility')
output_string_score += s1
 

s1 = similarity('lines_merged\\lines_merged_msavae_samples_ORIGINAL.txt', msa_outputs_samples, msa_reference, 'Samples for basic MSA-VAE model')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol0_samples_ORIGINAL.txt', msa_outputs_sol0_samples, msa_reference,'Samples for MSA-VAE model trained on low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol1_samples_ORIGINAL.txt', msa_outputs_sol1_samples, msa_reference,'Samples for MSA-VAE model trained on mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_sol2_samples_ORIGINAL.txt', msa_outputs_sol2_samples, msa_reference,'Samples for MSA-VAE model trained on high solubility')
output_string_score += s1
 

s1 = similarity('lines_merged\\lines_merged_arvae_variants_ORIGINAL.txt', ar_outputs_variants, ar_reference, 'Variants for basic AR-VAE model')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol0_variants_ORIGINAL.txt', ar_outputs_sol0_variants, ar_reference, 'Variants for AR-VAE model trained on low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol1_variants_ORIGINAL.txt', ar_outputs_sol1_variants, ar_reference, 'Variants for AR-VAE model trained on mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol2_variants_ORIGINAL.txt', ar_outputs_sol2_variants, ar_reference, 'Variants for AR-VAE model trained on high solubility')
output_string_score += s1
 

s1 = similarity('lines_merged\\lines_merged_arvae_samples_ORIGINAL.txt', ar_outputs_samples, ar_reference, 'Samples for basic AR-VAE model')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol0_samples_ORIGINAL.txt', ar_outputs_sol0_samples, ar_reference,' Samples for AR-VAE model trained on low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol1_samples_ORIGINAL.txt', ar_outputs_sol1_samples, ar_reference, 'Samples for AR-VAE model trained on mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_sol2_samples_ORIGINAL.txt', ar_outputs_sol2_samples, ar_reference, 'Samples for AR-VAE model trained on high solubility')
output_string_score += s1
 

s1 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol0_variants_ORIGINAL.txt', msa_outputs_cond_sol0, msa_reference_cond,'Variants for conditional MSA-VAE model with low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol1_variants_ORIGINAL.txt', msa_outputs_cond_sol1, msa_reference_cond,'Variants for conditional MSA-VAE model with mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_msavae_with_conditions_sol2_variants_ORIGINAL.txt', msa_outputs_cond_sol2, msa_reference_cond,'Variants for conditional MSA-VAE model with high solubility')
output_string_score += s1
 

s1 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol0_variants_ORIGINAL.txt', ar_outputs_cond_sol0, ar_reference_cond,'Variants for conditional AR-VAE model with low solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol1_variants_ORIGINAL.txt', ar_outputs_cond_sol1, ar_reference_cond,'Variants for conditional AR-VAE model with mid solubility')
output_string_score += s1
 
s1 = similarity('lines_merged\\lines_merged_arvae_with_conditions_sol2_variants_ORIGINAL.txt', ar_outputs_cond_sol2, ar_reference_cond,'Variants for conditional AR-VAE model with high solubility')
output_string_score += s1
 
   
# s1 = similarity('lines_merged\\lines_merged_PF00296_full.txt', 'Profile of luciferase family')
# output_string_score += s1
#  
s1 = similarity('training_validation\\luxafilt_llmsa_val.fa', ar_outputs, ar_reference, 'Validation data (AR)', '_ar')
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_val.fa', msa_outputs, msa_reference, 'Validation data (MSA)', '_msa')
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_val.fa', ar_outputs_cond, ar_reference_cond, 'Validation data (AR conditional)', '_ar_cond')
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_val.fa', msa_outputs_cond, msa_reference_cond, 'Validation data (MSA conditional)', '_msa_cond')
output_string_score += s1
 
s1 = similarity('training_validation\\luxafilt_llmsa_train.fa', ar_outputs, ar_reference, 'Training data (AR)', '_ar')
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_train.fa', msa_outputs, msa_reference, 'Training data (MSA)', '_msa') 
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_train.fa', ar_outputs_cond, ar_reference_cond, 'Training data (AR conditional)', '_ar_cond')
output_string_score += s1
s1 = similarity('training_validation\\luxafilt_llmsa_train.fa', msa_outputs_cond, msa_reference_cond, 'Training data (MSA conditional)', '_msa_cond')
output_string_score += s1

s1 = similarity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], ar_outputs, ar_reference, 'Training and validation data (AR)', '_ar')
output_string_score += s1
s1 = similarity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], msa_outputs, msa_reference, 'Training and validation data (MSA)', '_msa')
output_string_score += s1
s1 = similarity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], ar_outputs_cond, ar_reference_cond, 'Training and validation data (AR conditional)', '_ar_cond')
output_string_score += s1
s1 = similarity_multiple(['training_validation\\luxafilt_llmsa_train.fa', 'training_validation\\luxafilt_llmsa_val.fa'], msa_outputs_cond, msa_reference_cond, 'Training and validation data (MSA conditional)', '_msa_cond')
output_string_score += s1
 

output_string_score = output_string_score.replace(".", ",")
file_csv = open("..\\results\\tables\\cosine.csv", "w")
file_csv.write(output_string_score)
file_csv.close()