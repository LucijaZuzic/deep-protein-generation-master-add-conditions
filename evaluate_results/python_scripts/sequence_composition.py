# -*- coding: utf-8 -*-
"""
Created on Sun May 22 12:57:42 2022

@author: Lucija
"""

import matplotlib.pyplot as plt 
import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
positive = ['R','K','H']
negative = ['D','E']

all_mino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

hydro = {
    "I": -1.12,
    "L": -1.25,
    "F": -1.71,
    "V": -0.46,
    "M": -0.67,
    "P": 0.14,
    "W": -2.09,
    "H": 0.11,
    "T": 0.25,
    "Q": 0.77,
    "C": -0.02,
    "Y": -0.71,
    "A": 0.5,
    "S": 0.46,
    "N": 0.85,
    "R": 1.81,
    "G": 1.15,
    "E": 3.63,
    "K": 2.8,
    "D": 3.64,
}

def sequence_composition(filename, title):
    try:
        file = open("..\\sequences\\" + filename, 'r') 
    except:
        return 
    lines = file.readlines()
    lines = [line.strip().replace('-', '') for line in lines]
    new_charged = [] 
    new_hydro = [] 
    new_composition = {} 
    for amino_acid in all_mino_acids:
        new_composition[amino_acid] = 0
    for sequence_index in range(1, len(lines), 2): 
        line = lines[sequence_index]
        net_charge = 0
        net_hydro = 0
        for amino_acid in all_mino_acids:
            amino_count = line.count(amino_acid)
            new_composition[amino_acid] += line.count(amino_acid)
            net_hydro = hydro[amino_acid] * line.count(amino_acid)
            if amino_acid in positive:
                net_charge += amino_count
            elif amino_acid in negative:
                net_charge -= amino_count
        new_charged.append(IP(lines[sequence_index]).charge_at_pH(7.0)) 
        min_log_P = len(line) * min(hydro.values())
        max_log_P = len(line) * max(hydro.values())
        log_P_range = max_log_P - min_log_P
        new_hydro.append((net_hydro - min_log_P) / log_P_range)
    print(title) 
    print("charge", np.min(new_charged), np.max(new_charged), np.mean(new_charged), np.std(new_charged))
    plt.hist(new_charged)
    plt.xlabel("Charge")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\sequence_charge\\" + filename.replace('.txt', '_sequence_charge.png').replace('.fa', '_sequence_charge.png').replace('right_padded\\', '').replace('training_validation\\', '').replace('lines_merged\\lines_merged_', ''), bbox_inches='tight')
    plt.close()
    print("hydro", np.min(new_hydro), np.max(new_hydro), np.mean(new_hydro), np.std(new_hydro))
    plt.hist(new_hydro)
    plt.xlabel("log P")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\sequence_hydrophobicity\\" + filename.replace('.txt', '_sequence_hydrophobicity.png').replace('.fa', '_sequence_hydrophobicity.png').replace('right_padded\\', '').replace('training_validation\\', '').replace('lines_merged\\lines_merged_', ''), bbox_inches='tight')
    plt.close()
    number_amino_all = 0
    for amino_acid in all_mino_acids:
        number_amino_all += new_composition[amino_acid] 
    for amino_acid in all_mino_acids:
        new_composition[amino_acid] = new_composition[amino_acid] / number_amino_all
    x_axis = [i for i in range(len(all_mino_acids))]
    print("composition", new_composition)
    plt.bar(x_axis, new_composition.values())
    plt.xticks(x_axis, all_mino_acids)
    plt.xlabel("Aminacid")
    plt.ylabel("Percentage of occurences in sequence set")
    plt.title(title)
    plt.savefig("..\\results\\sequence_composition\\" + filename.replace('.txt', '_sequence_composition.png').replace('.fa', '_sequence_composition.png').replace('right_padded\\', '').replace('training_validation\\', '').replace('lines_merged\\lines_merged_', ''), bbox_inches='tight')
    plt.close()
    file.close()
    retval = ["","",""]
    retval[0] = title + ";" + str(np.round(np.min(new_charged),3)) + ";" + str(np.round(np.max(new_charged),3)) + ";" + str(np.round(np.mean(new_charged),3)) + ";" + str(np.round(np.std(new_charged),3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(new_hydro),3)) + ";" + str(np.round(np.max(new_hydro),3)) + ";" + str(np.round(np.mean(new_hydro),3)) + ";" + str(np.round(np.std(new_hydro),3))  + "\n"           
    retval[2] = title
    for amino_acid in all_mino_acids:
        retval[2] += ";" + str(np.round(new_composition[amino_acid] * 100,3))
    retval[2] += "\n"
    return retval
    
def sequence_composition_multiple(filenames, title):
    lines = []
    all_names = ""
    for filename in filenames:
        try:
            file = open("..\\sequences\\" + filename, 'r')
        except:
            return
        lines_part = file.readlines()
        all_names += filename.replace('.txt', '_').replace('.fa', '_').replace('right_padded\\', '').replace('training_validation\\', '').replace('lines_merged\\lines_merged_', '')
        lines_part = [line.strip().replace('-', '') for line in lines_part]
        lines += lines_part 
        file.close()
    new_charged = [] 
    new_hydro = [] 
    new_composition = {}
    for amino_acid in all_mino_acids:
        new_composition[amino_acid] = 0
    for sequence_index in range(1, len(lines), 2): 
        line = lines[sequence_index]
        net_charge = 0
        net_hydro = 0
        for amino_acid in all_mino_acids:
            amino_count = line.count(amino_acid)
            new_composition[amino_acid] += line.count(amino_acid)
            net_hydro = hydro[amino_acid] * line.count(amino_acid)
            if amino_acid in positive:
                net_charge += amino_count
            elif amino_acid in negative:
                net_charge -= amino_count
        new_charged.append(IP(lines[sequence_index]).charge_at_pH(7.0))
        min_log_P = len(line) * min(hydro.values())
        max_log_P = len(line) * max(hydro.values())
        log_P_range = max_log_P - min_log_P
        new_hydro.append((net_hydro - min_log_P) / log_P_range)
    print(title) 
    print("charge", np.min(new_charged), np.max(new_charged), np.mean(new_charged), np.std(new_charged))
    plt.hist(new_charged)
    plt.xlabel("Charge")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\sequence_charge\\" + all_names + "sequence_charge.png", bbox_inches='tight')
    plt.close()
    print("hydro", np.min(new_hydro), np.max(new_hydro), np.mean(new_hydro), np.std(new_hydro))
    plt.hist(new_hydro)
    plt.xlabel("log P")
    plt.ylabel("Number of sequences")
    plt.title(title)
    plt.savefig("..\\results\\sequence_hydrophobicity\\" + all_names + "sequence_hydrophobicity.png", bbox_inches='tight')
    plt.close()
    number_amino_all = 0
    for amino_acid in all_mino_acids:
        number_amino_all += new_composition[amino_acid] 
    for amino_acid in all_mino_acids:
        new_composition[amino_acid] = new_composition[amino_acid] / number_amino_all
    x_axis = [i for i in range(len(all_mino_acids))]
    print("composition", new_composition)
    plt.bar(x_axis, new_composition.values())
    plt.xticks(x_axis, all_mino_acids)
    plt.xlabel("Aminoacid")
    plt.ylabel("Percentage of occurences in sequence set")
    plt.title(title)
    plt.savefig("..\\results\\sequence_composition\\" + all_names + "sequence_composition.png", bbox_inches='tight')
    plt.close()
    file.close()
    retval = ["","",""]
    retval[0] = title + ";" + str(np.round(np.min(new_charged),3)) + ";" + str(np.round(np.max(new_charged),3)) + ";" + str(np.round(np.mean(new_charged),3)) + ";" + str(np.round(np.std(new_charged),3))  + "\n"       
    retval[1] = title + ";" + str(np.round(np.min(new_hydro),3)) + ";" + str(np.round(np.max(new_hydro),3)) + ";" + str(np.round(np.mean(new_hydro),3)) + ";" + str(np.round(np.std(new_hydro),3))  + "\n"       
    retval[2] = title
    for amino_acid in all_mino_acids:
        retval[2] += ";" + str(np.round(new_composition[amino_acid] * 100,3))
    retval[2] += "\n"
    return retval
  
output_string_charged = "Sequence set;Minimal charge;Maximal charge;Average charge;Standard deviation\n"
output_string_hydro = "Sequence set;Minimal log P;Maximal log P;Average log P;Standard deviation\n"
output_string_composition = ""
for amino_acid in all_mino_acids:
    output_string_composition += ";" + amino_acid
output_string_composition += "\n"

s1, s2, s3 = sequence_composition('right_padded\\msavae_variants.fa','Variants for basic MSA-VAE model')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol0_variants.fa','Variants for MSA-VAE model trained on low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol1_variants.fa','Variants for MSA-VAE model trained on mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol2_variants.fa','Variants for MSA-VAE model trained on high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

s1, s2, s3 = sequence_composition('right_padded\\msavae_samples.fa','Samples for basic MSA-VAE model')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol0_samples.fa','Samples for MSA-VAE model trained on low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol1_samples.fa','Samples for MSA-VAE model trained on mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_sol2_samples.fa','Samples for MSA-VAE model trained on high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
            
s1, s2, s3 = sequence_composition('right_padded\\arvae_variants.fa','Variants for basic AR-VAE model')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol0_variants.fa','Variants for AR-VAE model trained on low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol1_variants.fa','Variants for AR-VAE model trained on mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol2_variants.fa','Variants for AR-VAE model trained on high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

s1, s2, s3 = sequence_composition('right_padded\\arvae_samples.fa','Samples for basic AR-VAE model')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol0_samples.fa','Samples for AR-VAE model trained on low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol1_samples.fa','Samples for AR-VAE model trained on mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_sol2_samples.fa','Samples for AR-VAE model trained on high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

s1, s2, s3 = sequence_composition('right_padded\\msavae_with_conditions_sol0_variants.fa','Variants for conditional MSA-VAE model with low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_with_conditions_sol1_variants.fa','Variants for conditional MSA-VAE model with mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\msavae_with_conditions_sol2_variants.fa','Variants for conditional MSA-VAE model with high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

s1, s2, s3 = sequence_composition('right_padded\\arvae_with_conditions_sol0_variants.fa','Variants for conditional AR-VAE model with low solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_with_conditions_sol1_variants.fa','Variants for conditional AR-VAE model with mid solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('right_padded\\arvae_with_conditions_sol2_variants.fa','Variants for conditional AR-VAE model with high solubility')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

s1, s2, s3 = sequence_composition('lines_merged\\lines_merged_PF00296_full.txt', 'Profile of luciferase family')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('training_validation\\ll_val.fa', 'Validation data')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition('training_validation\\ll_train.fa', 'Training data')  
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3
s1, s2, s3 = sequence_composition_multiple(['training_validation\\ll_train.fa', 'training_validation\\ll_val.fa'], 'Training and validation data')
output_string_charged += s1
output_string_hydro += s2
output_string_composition += s3

output_string_hydro = output_string_hydro.replace(".", ",")
file_csv = open("..\\results\\tables\\sequence_hydrophobicity.csv", "w")
file_csv.write(output_string_hydro)
file_csv.close()

output_string_charged = output_string_charged.replace(".", ",")
file_csv = open("..\\results\\tables\\sequence_charge.csv", "w")
file_csv.write(output_string_charged)
file_csv.close()

output_string_composition = output_string_composition.replace(".", ",")
file_csv = open("..\\results\\tables\\sequence_composition.csv", "w")
file_csv.write(output_string_composition)
file_csv.close()