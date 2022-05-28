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
    plt.xlabel("Skalirana predviđena topljivost")
    plt.ylabel("Broj sekvenci")
    plt.axvline(0.446, color='red', linestyle='dotted', label="Referentna topljivost") 
    plt.axvline(0.447, color='green', linestyle='dashed', label="Donja granica visoke topljivosti")
    plt.axvline(0.34, color='orange', linestyle='dashed', label="Donja granica srednje topljivosti") 
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
    plt.xlabel("Skalirana predviđena topljivost")
    plt.ylabel("Broj sekvenci")
    plt.axvline(0.446, color='red', linestyle='dotted', label="Referentna topljivost") 
    plt.axvline(0.447, color='green', linestyle='dashed', label="Donja granica visoke topljivosti")
    plt.axvline(0.34, color='orange', linestyle='dashed', label="Donja granica srednje topljivosti") 
    plt.legend(bbox_to_anchor=(0, -0.4), loc="lower left")
    plt.title(title)
    plt.savefig("..\\results\\read_prediction\\" + all_names + 'multiple.png', bbox_inches='tight') 
    plt.close()
    return title + ";" + str(np.round(np.min(sol),3)) + ";" + str(np.round(np.max(sol),3)) + ";" + str(np.round(np.mean(sol),3)) + ";" + str(np.round(np.std(sol),3))  + "\n"  

output_string = "Skup sekvenci;Minimalna topljivost;Maksimalna topljivost;Prosjek topljivosti;Standardna devijacija\n"

output_string += read_prediction('seq_prediction_msavae_variants.txt','Varijante za osnovni MSA-VAE model')
output_string += read_prediction('seq_prediction_msavae_sol0_variants.txt','Varijante za MSA-VAE model niske topljivosti')
output_string += read_prediction('seq_prediction_msavae_sol1_variants.txt','Varijante za MSA-VAE model srednje topljivosti')
output_string += read_prediction('seq_prediction_msavae_sol2_variants.txt','Varijante za MSA-VAE model visoke topljivosti')
output_string += read_prediction('seq_prediction_msavae_samples.txt','Uzorci za osnovni MSA-VAE model')
output_string += read_prediction('seq_prediction_msavae_sol0_samples.txt','Uzorci za MSA-VAE model niske topljivosti')
output_string += read_prediction('seq_prediction_msavae_sol1_samples.txt','Uzorci za MSA-VAE model srednje topljivosti')
output_string += read_prediction('seq_prediction_msavae_sol2_samples.txt','Uzorci za MSA-VAE model visoke topljivosti')
            
output_string += read_prediction('seq_prediction_arvae_variants.txt','Varijante za osnovni AR-VAE model')
output_string += read_prediction('seq_prediction_arvae_sol0_variants.txt','Varijante za AR-VAE model niske topljivosti')
output_string += read_prediction('seq_prediction_arvae_sol1_variants.txt','Varijante za AR-VAE model srednje topljivosti')
output_string += read_prediction('seq_prediction_arvae_sol2_variants.txt','Varijante za AR-VAE model visoke topljivosti')
output_string += read_prediction('seq_prediction_arvae_samples.txt','Uzorci za osnovni AR-VAE model')
output_string += read_prediction('seq_prediction_arvae_sol0_samples.txt','Uzorci za AR-VAE model niske topljivosti')
output_string += read_prediction('seq_prediction_arvae_sol1_samples.txt','Uzorci za AR-VAE model srednje topljivosti')
output_string += read_prediction('seq_prediction_arvae_sol2_samples.txt','Uzorci za AR-VAE model visoke topljivosti')

output_string += read_prediction('seq_prediction_msavae_with_conditions_sol0_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom niske topljivosti')
output_string += read_prediction('seq_prediction_msavae_with_conditions_sol1_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom srednje topljivosti')
output_string += read_prediction('seq_prediction_msavae_with_conditions_sol2_variants.txt','Varijante za uvjetni MSA-VAE model s uvjetom visoke topljivosti')

output_string += read_prediction('seq_prediction_arvae_with_conditions_sol0_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom niske topljivosti')
output_string += read_prediction('seq_prediction_arvae_with_conditions_sol1_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom srednje topljivosti')
output_string += read_prediction('seq_prediction_arvae_with_conditions_sol2_variants.txt','Varijante za uvjetni AR-VAE model s uvjetom visoke topljivosti')

output_string += read_prediction('seq_prediction_PF00296_full.txt', 'Profil obitelji luciferaza')
output_string += read_prediction('seq_prediction_ll_val.txt', 'Skup podataka za validaciju')
output_string += read_prediction('seq_prediction_ll_train.txt', 'Skup podataka za treniranje')

output_string += read_prediction_multiple(['seq_prediction_ll_train.txt', 'seq_prediction_ll_val.txt'], 'Skup podataka za treniranje i validaciju')

output_string = output_string.replace(".", ",")
file_csv = open("..\\results\\tables\\read_prediction.csv", "w")
file_csv.write(output_string)
file_csv.close()