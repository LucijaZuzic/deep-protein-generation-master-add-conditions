# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:54:40 2022

@author: Lucija
"""
import numpy as np 
 
def true_prediction(filename, title, limit1, limit2, bin_number):
    try:
        file = open("../sequences/seq_prediction/" + filename, 'r')
    except:
        return
    lines = file.readlines()
    bins = [0,0,0]
    for line in lines: 
        processed_line = line.replace(' ', '').replace('\n', '').split(',')
        if 'SEQUENCEPREDICTIONS' == processed_line[0]:  
            if float(processed_line[3]) < limit1:
                bins[0] += 1
                continue
            if float(processed_line[3]) < limit2:
                bins[1] += 1
                continue
            bins[2] += 1 
    sum_predictions = np.sum(bins)
    bins = [one_bin / sum_predictions for one_bin in bins]
    true_predictions = bins[bin_number]
    print(title, bins, true_predictions)       
    file.close()
    return title + ";" + str(np.round(bins[0] * 100, 3)) + ";" + str(np.round(bins[1] * 100, 3)) + ";" + str(np.round(bins[2] * 100, 3)) + "\n"
    
def class_division(filename, title, limit1, limit2):
    try:
        file = open("../sequences/seq_prediction/" + filename, 'r')
    except:
        return
    lines = file.readlines()
    bins = [0,0,0]
    for line in lines: 
        processed_line = line.replace(' ', '').replace('\n', '').split(',')
        if 'SEQUENCEPREDICTIONS' == processed_line[0]:  
            if float(processed_line[3]) < limit1:
                bins[0] += 1
                continue
            if float(processed_line[3]) < limit2:
                bins[1] += 1
                continue
            bins[2] += 1 
    sum_predictions = np.sum(bins)
    bins = [one_bin / sum_predictions for one_bin in bins]
    print(title, bins)       
    file.close()
    return title + ";" + str(np.round(bins[0] * 100, 3)) + ";" + str(np.round(bins[1] * 100, 3)) + ";" + str(np.round(bins[2] * 100, 3)) + "\n"
    
def class_division_multiple(filenames, title, limit1, limit2):
    lines = []
    for filename in filenames:
        try:
            file = open("../sequences/seq_prediction/" + filename, 'r')
        except:
            return
        lines += file.readlines() 
        file.close() 
    bins = [0,0,0]
    for line in lines: 
        processed_line = line.replace(' ', '').replace('\n', '').split(',')
        if 'SEQUENCEPREDICTIONS' == processed_line[0]:  
            if float(processed_line[3]) < limit1:
                bins[0] += 1
                continue
            if float(processed_line[3]) < limit2:
                bins[1] += 1
                continue
            bins[2] += 1 
    sum_predictions = np.sum(bins)
    bins = [one_bin / sum_predictions for one_bin in bins]
    print(title, bins)       
    file.close()
    return title + ";" + str(np.round(bins[0] * 100, 3)) + ";" + str(np.round(bins[1] * 100, 3)) + ";" + str(np.round(bins[2] * 100, 3)) + "\n"
    
output_string =  "Sequence set;Percentage of low solubility samples;Percentage of mid solubility samples;Percentage of high solubility samples\n"
output_string += class_division('seq_prediction_msavae_variants.txt','Variants for basic MSA-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_msavae_sol0_variants.txt','Variants for MSA-VAE model trained on low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_sol1_variants.txt','Variants for MSA-VAE model trained on mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_sol2_variants.txt','Variants for MSA-VAE model trained on high solubility', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_msavae_samples.txt','Samples for MSA-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_msavae_sol0_samples.txt','Samples for MSA-VAE model trained on low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_sol1_samples.txt','Samples for MSA-VAE model trained on mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_sol2_samples.txt','Samples for MSA-VAE model trained on high solubility', 0.34, 0.447, 2)
            
output_string += class_division('seq_prediction_arvae_variants.txt','Variants for for AR-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_arvae_sol0_variants.txt','Variants for AR-VAE model trained on low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_sol1_variants.txt','Variants for AR-VAE model trained on mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_sol2_variants.txt','Variants for AR-VAE model trained on high solubility', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_arvae_samples.txt','Samples for AR-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_arvae_sol0_samples.txt','Samples for AR-VAE model trained on low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_sol1_samples.txt','Samples for AR-VAE model trained on mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_sol2_samples.txt','Samples for AR-VAE model trained on high solubility', 0.34, 0.447, 2)

output_string += true_prediction('seq_prediction_msavae_with_conditions_sol0_variants.txt','Variants for conditional MSA-VAE model with low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_with_conditions_sol1_variants.txt','Variants for conditional MSA-VAE model with mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_with_conditions_sol2_variants.txt','Variants for conditional MSA-VAE model with high solubility', 0.34, 0.447, 2)

output_string += true_prediction('seq_prediction_arvae_with_conditions_sol0_variants.txt','Variants for conditional AR-VAE model with low solubility', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_with_conditions_sol1_variants.txt','Variants for conditional AR-VAE model with mid solubility', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_with_conditions_sol2_variants.txt','Variants for conditional AR-VAE model with high solubility', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_PF00296_full.txt', 'Profile of luciferase family', 0.34, 0.447)
output_string += class_division('seq_prediction_ll_val.txt', 'Validation data', 0.34, 0.447)
output_string += class_division('seq_prediction_ll_train.txt', 'Training data', 0.34, 0.447)

output_string += class_division_multiple(['seq_prediction_ll_train.txt', 'seq_prediction_ll_val.txt'], 'Training and validation data', 0.34, 0.447)

output_string = output_string.replace(".", ",")
file_csv = open("../results/tables/true_prediction.csv", "w")
file_csv.write(output_string)
file_csv.close()