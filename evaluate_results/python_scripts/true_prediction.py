# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:54:40 2022

@author: Lucija
"""
import numpy as np 
 
def true_prediction(filename, title, limit1, limit2, bin_number):
    try:
        file = open("..\\sequences\\seq_prediction\\" + filename, 'r')
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
        file = open("..\\sequences\\seq_prediction\\" + filename, 'r')
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
            file = open("..\\sequences\\seq_prediction\\" + filename, 'r')
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
    
output_string =  "Skup sekvenci;Udio uzoraka niske topljivosti (postotak);Udio uzoraka srednje topljivosti (postotak);Udio uzoraka visoke topljivosti (postotak)\n"
output_string += class_division('seq_prediction_msavae_variants.txt','Varijante za osnovni MSA-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_msavae_sol0_variants.txt','Varijante za MSA-VAE model niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_sol1_variants.txt','Varijante za MSA-VAE model srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_sol2_variants.txt','Varijante za MSA-VAE model visoke topljivosti', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_msavae_samples.txt','Uzorci za osnovni MSA-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_msavae_sol0_samples.txt','Uzorci za MSA-VAE model niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_sol1_samples.txt','Uzorci za MSA-VAE model srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_sol2_samples.txt','Uzorci za MSA-VAE model visoke topljivosti', 0.34, 0.447, 2)
            
output_string += class_division('seq_prediction_arvae_variants.txt','Varijante za osnovni AR-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_arvae_sol0_variants.txt','Varijante za AR-VAE model niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_sol1_variants.txt','Varijante za AR-VAE model srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_sol2_variants.txt','Varijante za AR-VAE model visoke topljivosti', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_arvae_samples.txt','Uzorci za osnovni AR-VAE model', 0.34, 0.447)
output_string += true_prediction('seq_prediction_arvae_sol0_samples.txt','Uzorci za AR-VAE model niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_sol1_samples.txt','Uzorci za AR-VAE model srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_sol2_samples.txt','Uzorci za AR-VAE model visoke topljivosti', 0.34, 0.447, 2)

output_string += true_prediction('seq_prediction_msavae_with_conditions_sol0_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_msavae_with_conditions_sol1_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_msavae_with_conditions_sol2_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom visoke topljivosti', 0.34, 0.447, 2)

output_string += true_prediction('seq_prediction_arvae_with_conditions_sol0_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom niske topljivosti', 0.34, 0.447, 0)
output_string += true_prediction('seq_prediction_arvae_with_conditions_sol1_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom srednje topljivosti', 0.34, 0.447, 1)
output_string += true_prediction('seq_prediction_arvae_with_conditions_sol2_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom visoke topljivosti', 0.34, 0.447, 2)

output_string += class_division('seq_prediction_PF00296_full.txt', 'Profil obitelji luciferaza', 0.34, 0.447)
output_string += class_division('seq_prediction_ll_val.txt', 'Skup podataka za validaciju', 0.34, 0.447)
output_string += class_division('seq_prediction_ll_train.txt', 'Skup podataka za treniranje', 0.34, 0.447)

output_string += class_division_multiple(['seq_prediction_ll_train.txt', 'seq_prediction_ll_val.txt'], 'Skup podataka za treniranje i validaciju', 0.34, 0.447)

output_string = output_string.replace(".", ",")
file_csv = open("..\\results\\tables\\true_prediction.csv", "w")
file_csv.write(output_string)
file_csv.close()