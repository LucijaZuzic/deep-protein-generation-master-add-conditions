# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 18:54:27 2022

@author: Lucija
"""
import numpy as np 
import matplotlib.pyplot as plt
BLOSSUM_MATRIX = []
characters = []


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
 
def find_letter(letter):
    for i in range(len(characters)):
        if characters[i] == letter:
            return i
    return -1
    

aaletters = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def letter_score(scores, name):
    for index_a in range(1, len(aaletters)):
        for index_b in range(index_a + 1, len(aaletters)): 
            cosine = np.dot(scores[index_a],scores[index_b])/(np.linalg.norm(scores[index_a])*np.linalg.norm(scores[index_b]))
            blossum = BLOSSUM_MATRIX[find_letter(aaletters[index_a])][find_letter(aaletters[index_b])]
            plt.scatter(cosine, blossum, color="blue")
    plt.ylabel("BLOSSUM62 score")
    plt.xlabel("Cosine similarity")  
    plt.title("Amino acid pairs")  
    plt.savefig("../results/cosine/" + name + "_cosine_letter.png", bbox_inches='tight') 
    plt.close() 
              
msa_outputs_aaletters = np.load('../../output/generated_PCA/msa_outputs_aaletters.npy')[0]
letter_score(msa_outputs_aaletters, 'msa')
msa_outputs_aaletters_sol0 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol0.npy')[0]
letter_score(msa_outputs_aaletters_sol0, 'msa_sol0')
msa_outputs_aaletters_sol1 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol1.npy')[0]
letter_score(msa_outputs_aaletters_sol1, 'msa_sol1')
msa_outputs_aaletters_sol2 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol2.npy')[0]
letter_score(msa_outputs_aaletters_sol2, 'msa_sol2')
msa_outputs_aaletters_cond_sol0 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol0.npy')[0]
letter_score(msa_outputs_aaletters_cond_sol0, 'msa_cond_sol0')
msa_outputs_aaletters_cond_sol1 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol1.npy')[0]
letter_score(msa_outputs_aaletters_cond_sol1, 'msa_cond_sol1')
msa_outputs_aaletters_cond_sol2 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol2.npy')[0]
letter_score(msa_outputs_aaletters_cond_sol2, 'msa_cond_sol2')

ar_outputs_aaletters = np.load('../../output/generated_PCA/ar_outputs_aaletters.npy')[0]
letter_score(ar_outputs_aaletters, 'ar')
ar_outputs_aaletters_sol0 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol0.npy')[0]
letter_score(ar_outputs_aaletters_sol0, 'ar_sol0')
ar_outputs_aaletters_sol1 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol1.npy')[0]
letter_score(ar_outputs_aaletters_sol1, 'ar_sol1')
ar_outputs_aaletters_sol2 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol2.npy')[0]
letter_score(ar_outputs_aaletters_sol2, 'ar_sol2')
ar_outputs_aaletters_cond_sol0 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol0.npy')[0]
letter_score(ar_outputs_aaletters_cond_sol0, 'ar_cond_sol0')
ar_outputs_aaletters_cond_sol1 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol1.npy')[0]
letter_score(ar_outputs_aaletters_cond_sol1, 'ar_cond_sol1')
ar_outputs_aaletters_cond_sol2 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol2.npy')[0]
letter_score(ar_outputs_aaletters_cond_sol2, 'ar_cond_sol2')
