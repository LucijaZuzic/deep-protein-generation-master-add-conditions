# -*- coding: utf-8 -*-
"""
Created on Sun May 22 12:57:42 2022

@author: Lucija
"""
import matplotlib.pyplot as plt
import numpy as np

def read_prediction(filename, title):
    try:
        file = open("..\\sequences\\seq_prediction\\" + filename, 'r')
    except:
        return
    lines = file.readlines()
    sol = []
    for line in lines: 
        processed_line = line.replace(' ', '').replace('\n', '').split(',')
        if 'SEQUENCEPREDICTIONS' == processed_line[0]: 
            sol.append(float(processed_line[3]))
    print(title, np.min(sol), np.max(sol), np.mean(sol), np.std(sol))
    plt.hist(sol)
    plt.xlabel("Scaled predicted solubility") 
    plt.ylabel("Number of sequences")
    plt.axvline(0.446, color='red', linestyle='dotted', label="Reference solubility")
    plt.axvline(0.447, color='green', linestyle='dashed', label="Lower bound of high solubility")
    plt.axvline(0.34, color='orange', linestyle='dashed', label="Lower bound of mid solubility") 
    plt.legend(bbox_to_anchor=(0, -0.4), loc="lower left")
    plt.title(title)
    plt.savefig("..\\results\\read_prediction\\" + filename.replace('txt', 'png'), bbox_inches='tight')
    plt.close()
    file.close()
    return title + ";" + str(np.round(np.min(sol),3)) + ";" + str(np.round(np.max(sol),3)) + ";" + str(np.round(np.mean(sol),3)) + ";" + str(np.round(np.std(sol),3))  + "\n"       
    
def read_prediction_multiple(filenames, title):
    lines = []
    all_names = ""
    for filename in filenames:
        try:
            file = open("..\\sequences\\seq_prediction\\" + filename, 'r')
        except:
            return
        lines += file.readlines()
        all_names += filename.replace('.txt', '_')
        file.close()
    sol = []
    for line in lines: 
        processed_line = line.replace(' ', '').replace('\n', '').split(',')
        if 'SEQUENCEPREDICTIONS' == processed_line[0]: 
            sol.append(float(processed_line[3]))
    print(title, np.min(sol), np.max(sol), np.mean(sol), np.std(sol))
    plt.hist(sol)
    plt.xlabel("Scaled predicted solubility")
    plt.ylabel("Number of sequences")
    plt.axvline(0.446, color='red', linestyle='dotted', label="Reference solubility") 
    plt.axvline(0.447, color='green', linestyle='dashed', label="Lower bound of high solubility")
    plt.axvline(0.34, color='orange', linestyle='dashed', label="Lower bound of mid solubility") 
    plt.legend(bbox_to_anchor=(0, -0.4), loc="lower left")
    plt.title(title)
    plt.savefig("..\\results\\read_prediction\\" + all_names + 'multiple.png', bbox_inches='tight') 
    plt.close()
    return title + ";" + str(np.round(np.min(sol),3)) + ";" + str(np.round(np.max(sol),3)) + ";" + str(np.round(np.mean(sol),3)) + ";" + str(np.round(np.std(sol),3))  + "\n"  

output_string = "Sequence set;Minimal solubility;Maximal solubility;Average solubility;Standard deviation\n"

output_string += read_prediction('seq_prediction_msavae_variants.txt','Variants for basic MSA-VAE model')
output_string += read_prediction('seq_prediction_msavae_sol0_variants.txt','Variants for MSA-VAE model trained on low solubility')
output_string += read_prediction('seq_prediction_msavae_sol1_variants.txt','Variants for MSA-VAE model trained on mid solubility')
output_string += read_prediction('seq_prediction_msavae_sol2_variants.txt','Variants for MSA-VAE model trained on high solubility')
output_string += read_prediction('seq_prediction_msavae_samples.txt','Samples for basic MSA-VAE model')
output_string += read_prediction('seq_prediction_msavae_sol0_samples.txt','Samples for MSA-VAE model trained on low solubility')
output_string += read_prediction('seq_prediction_msavae_sol1_samples.txt','Samples for MSA-VAE model trained on mid solubility')
output_string += read_prediction('seq_prediction_msavae_sol2_samples.txt','Samples for MSA-VAE model trained on high solubility')
            
output_string += read_prediction('seq_prediction_arvae_variants.txt','Variants for basic AR-VAE model')
output_string += read_prediction('seq_prediction_arvae_sol0_variants.txt','Variants for AR-VAE model trained on low solubility')
output_string += read_prediction('seq_prediction_arvae_sol1_variants.txt','Variants for AR-VAE model trained on mid solubility')
output_string += read_prediction('seq_prediction_arvae_sol2_variants.txt','Variants fo rAR-VAE model trained on high solubility')
output_string += read_prediction('seq_prediction_arvae_samples.txt','Samples for basic AR-VAE model')
output_string += read_prediction('seq_prediction_arvae_sol0_samples.txt','Samples for AR-VAE model trained on low solubility')
output_string += read_prediction('seq_prediction_arvae_sol1_samples.txt','Samples for AR-VAE model trained on mid solubility')
output_string += read_prediction('seq_prediction_arvae_sol2_samples.txt','Samples for AR-VAE model trained on high solubility')

output_string += read_prediction('seq_prediction_msavae_with_conditions_sol0_variants.txt','Variants for conditional MSA-VAE model with low solubility')
output_string += read_prediction('seq_prediction_msavae_with_conditions_sol1_variants.txt','Variants for conditional MSA-VAE model with mid solubility')
output_string += read_prediction('seq_prediction_msavae_with_conditions_sol2_variants.txt','Variants for conditional MSA-VAE model with high solubility')

output_string += read_prediction('seq_prediction_arvae_with_conditions_sol0_variants.txt','Variants for conditional AR-VAE model with low solubility')
output_string += read_prediction('seq_prediction_arvae_with_conditions_sol1_variants.txt','Variants for conditional AR-VAE model with mid solubility')
output_string += read_prediction('seq_prediction_arvae_with_conditions_sol2_variants.txt','Variants for conditional AR-VAE model with high solubility')

output_string += read_prediction('seq_prediction_PF00296_full.txt', 'Profile of luciferase family') 
output_string += read_prediction('seq_prediction_ll_val.txt', 'Validation data')
output_string += read_prediction('seq_prediction_ll_train.txt', 'Training data')

output_string += read_prediction_multiple(['seq_prediction_ll_train.txt', 'seq_prediction_ll_val.txt'], 'Training and validation data')

output_string = output_string.replace(".", ",")
file_csv = open("..\\results\\tables\\read_prediction.csv", "w")
file_csv.write(output_string)
file_csv.close()