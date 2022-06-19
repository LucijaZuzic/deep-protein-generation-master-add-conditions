# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 12:30:11 2022

@author: Lucija
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

aromatic = ['H','Y','F','W']
negative = ['D','E']
positive = ['R','K']
hydrophobic = ['L','I','M','V','A']
polar = ['T','S','N','Q']
unique = ['C','P','G']
aaletters = aromatic + negative + positive + hydrophobic + polar + unique

def plot_letter(x, y, letter):
	if letter in aromatic:
		if letter == aromatic[0]:
			plt.scatter(x, y, color='blue', label = 'Aromatic')
		else:
			plt.scatter(x, y, color='blue')
	if letter in negative:
		if letter == negative[0]:
			plt.scatter(x, y, color='orange', label = 'Negative')
		else:
			plt.scatter(x, y, color='orange')
	if letter in positive:
		if letter == positive[0]:
			plt.scatter(x, y, color='green', label = 'Positive')
		else:
			plt.scatter(x, y, color='green') 
	if letter in hydrophobic:
		if letter == hydrophobic[0]:
			plt.scatter(x, y, color='red', label = 'Hydrophobic')
		else:
			plt.scatter(x, y, color='red')
	if letter in polar:
		if letter == polar[0]:
			plt.scatter(x, y, color='purple', label = 'Polar')
		else:
			plt.scatter(x, y, color='purple')
	if letter in unique:
		if letter == unique[0]:
			plt.scatter(x, y, color='brown', label = 'Unique')
		else:
			plt.scatter(x, y, color='brown')
            
pca = PCA(n_components=2)
 
outputs_aaletters = np.load('../../output/generated_PCA/msa_outputs_aaletters.npy')

outputs_aaletters_sol0 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol0.npy')
outputs_aaletters_sol1 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol1.npy')
outputs_aaletters_sol2 = np.load('../../output/generated_PCA/msa_outputs_aaletters_sol2.npy') 

outputs_aaletters_cond_sol0 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol0.npy')
outputs_aaletters_cond_sol1 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol1.npy')
outputs_aaletters_cond_sol2 = np.load('../../output/generated_PCA/msa_outputs_aaletters_cond_sol2.npy') 

outputs = np.load('../../output/generated_PCA/msa_outputs.npy')
outputs_samples = np.load('../../output/generated_PCA/msa_outputs_samples.npy')
outputs_variants = np.load('../../output/generated_PCA/msa_outputs_variants.npy')

outputs_sol0 = np.load('../../output/generated_PCA/msa_outputs_sol0.npy')
outputs_sol0_samples = np.load('../../output/generated_PCA/msa_outputs_sol0_samples.npy')
outputs_sol0_variants = np.load('../../output/generated_PCA/msa_outputs_sol0_variants.npy')

outputs_sol1 = np.load('../../output/generated_PCA/msa_outputs_sol1.npy')
outputs_sol1_samples = np.load('../../output/generated_PCA/msa_outputs_sol1_samples.npy')
outputs_sol1_variants = np.load('../../output/generated_PCA/msa_outputs_sol1_variants.npy')

outputs_sol2 = np.load('../../output/generated_PCA/msa_outputs_sol2.npy')
outputs_sol2_samples = np.load('../../output/generated_PCA/msa_outputs_sol2_samples.npy')
outputs_sol2_variants = np.load('../../output/generated_PCA/msa_outputs_sol2_variants.npy')

outputs_cond = np.load('../../output/generated_PCA/msa_outputs_cond.npy')
outputs_cond_sol0 = np.load('../../output/generated_PCA/msa_outputs_cond_sol0.npy')
outputs_cond_sol1 = np.load('../../output/generated_PCA/msa_outputs_cond_sol1.npy')
outputs_cond_sol2 = np.load('../../output/generated_PCA/msa_outputs_cond_sol2.npy')

pca_features1_aaletters = pca.fit_transform(outputs_aaletters[0]).reshape(2, -1)
pca_features2_aaletters = pca.fit_transform(outputs_aaletters[1]).reshape(2, -1)

pca_features1_aaletters_sol0 = pca.fit_transform(outputs_aaletters_sol0[0]).reshape(2, -1)
pca_features2_aaletters_sol0 = pca.fit_transform(outputs_aaletters_sol0[1]).reshape(2, -1)

pca_features1_aaletters_sol1 = pca.fit_transform(outputs_aaletters_sol1[0]).reshape(2, -1)
pca_features2_aaletters_sol1 = pca.fit_transform(outputs_aaletters_sol1[1]).reshape(2, -1)

pca_features1_aaletters_sol2 = pca.fit_transform(outputs_aaletters_sol2[0]).reshape(2, -1)
pca_features2_aaletters_sol2 = pca.fit_transform(outputs_aaletters_sol2[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol0 = pca.fit_transform(outputs_aaletters_cond_sol0[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol0 = pca.fit_transform(outputs_aaletters_cond_sol0[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol1 = pca.fit_transform(outputs_aaletters_cond_sol1[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol1 = pca.fit_transform(outputs_aaletters_cond_sol1[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol2 = pca.fit_transform(outputs_aaletters_cond_sol2[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol2 = pca.fit_transform(outputs_aaletters_cond_sol2[1]).reshape(2, -1)

pca_features1 = pca.fit_transform(outputs[0]).reshape(2, -1)
pca_features2 = pca.fit_transform(outputs[1]).reshape(2, -1)

pca_features1_samples = pca.fit_transform(outputs_samples[0]).reshape(2, -1)
pca_features2_samples = pca.fit_transform(outputs_samples[1]).reshape(2, -1)

pca_features1_variants = pca.fit_transform(outputs_variants[0]).reshape(2, -1)
pca_features2_variants = pca.fit_transform(outputs_variants[1]).reshape(2, -1)

pca_features1_sol0 = pca.fit_transform(outputs_sol0[0]).reshape(2, -1)
pca_features2_sol0 = pca.fit_transform(outputs_sol0[1]).reshape(2, -1)

pca_features1_sol0_samples = pca.fit_transform(outputs_sol0_samples[0]).reshape(2, -1)
pca_features2_sol0_samples = pca.fit_transform(outputs_sol0_samples[1]).reshape(2, -1)

pca_features1_sol0_variants = pca.fit_transform(outputs_sol0_variants[0]).reshape(2, -1)
pca_features2_sol0_variants = pca.fit_transform(outputs_sol0_variants[1]).reshape(2, -1)

pca_features1_sol1 = pca.fit_transform(outputs_sol1[0]).reshape(2, -1)
pca_features2_sol1 = pca.fit_transform(outputs_sol1[1]).reshape(2, -1)

pca_features1_sol1_samples = pca.fit_transform(outputs_sol1_samples[0]).reshape(2, -1)
pca_features2_sol1_samples = pca.fit_transform(outputs_sol1_samples[1]).reshape(2, -1)

pca_features1_sol1_variants = pca.fit_transform(outputs_sol1_variants[0]).reshape(2, -1)
pca_features2_sol1_variants = pca.fit_transform(outputs_sol1_variants[1]).reshape(2, -1)

pca_features1_sol2 = pca.fit_transform(outputs_sol2[0]).reshape(2, -1)
pca_features2_sol2 = pca.fit_transform(outputs_sol2[1]).reshape(2, -1)

pca_features1_sol2_samples = pca.fit_transform(outputs_sol2_samples[0]).reshape(2, -1)
pca_features2_sol2_samples = pca.fit_transform(outputs_sol2_samples[1]).reshape(2, -1)

pca_features1_sol2_variants = pca.fit_transform(outputs_sol2_variants[0]).reshape(2, -1)
pca_features2_sol2_variants = pca.fit_transform(outputs_sol2_variants[1]).reshape(2, -1)

pca_features1_cond = pca.fit_transform(outputs_cond[0]).reshape(2, -1)
pca_features2_cond = pca.fit_transform(outputs_cond[1]).reshape(2, -1)

pca_features1_sol0_cond = pca.fit_transform(outputs_cond_sol0[0]).reshape(2, -1)
pca_features2_sol0_cond = pca.fit_transform(outputs_cond_sol0[1]).reshape(2, -1)

pca_features1_sol1_cond = pca.fit_transform(outputs_cond_sol1[0]).reshape(2, -1)
pca_features2_sol1_cond = pca.fit_transform(outputs_cond_sol1[1]).reshape(2, -1)

pca_features1_sol2_cond = pca.fit_transform(outputs_cond_sol2[0]).reshape(2, -1)
pca_features2_sol2_cond = pca.fit_transform(outputs_cond_sol2[1]).reshape(2, -1)

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters[0][letter_index], pca_features1_aaletters[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters[0][letter_index], pca_features1_aaletters[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters[0][letter_index], pca_features2_aaletters[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters[0][letter_index], pca_features2_aaletters[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol0[0][letter_index], pca_features1_aaletters_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol0[0][letter_index], pca_features1_aaletters_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol0_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol0[0][letter_index], pca_features2_aaletters_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol0[0][letter_index], pca_features2_aaletters_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol0_PCA2.png', bbox_inches='tight') 
plt.close()
 
for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol1[0][letter_index], pca_features1_aaletters_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol1[0][letter_index], pca_features1_aaletters_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol1_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol1[0][letter_index], pca_features2_aaletters_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol1[0][letter_index], pca_features2_aaletters_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol1_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol2[0][letter_index], pca_features1_aaletters_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol2[0][letter_index], pca_features1_aaletters_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol2_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol2[0][letter_index], pca_features2_aaletters_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol2[0][letter_index], pca_features2_aaletters_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_sol2_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol0[0][letter_index], pca_features1_aaletters_cond_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol0[0][letter_index], pca_features1_aaletters_cond_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol0_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol0[0][letter_index], pca_features2_aaletters_cond_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol0[0][letter_index], pca_features2_aaletters_cond_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol0_PCA2.png', bbox_inches='tight') 
plt.close()
 
for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol1[0][letter_index], pca_features1_aaletters_cond_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol1[0][letter_index], pca_features1_aaletters_cond_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol1_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol1[0][letter_index], pca_features2_aaletters_cond_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol1[0][letter_index], pca_features2_aaletters_cond_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol1_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol2[0][letter_index], pca_features1_aaletters_cond_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol2[0][letter_index], pca_features1_aaletters_cond_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol2_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol2[0][letter_index], pca_features2_aaletters_cond_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol2[0][letter_index], pca_features2_aaletters_cond_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_aaletters_cond_sol2_PCA2.png', bbox_inches='tight') 
plt.close()
 
plt.scatter(pca_features1[0], pca_features1[1], label = 'all data')
plt.scatter(pca_features1_sol0[0], pca_features1_sol0[1], label = 'low solubility')
plt.scatter(pca_features1_sol1[0], pca_features1_sol1[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2[0], pca_features1_sol2[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for training and validation dataset split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_split_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'all data')
plt.scatter(pca_features2_sol0[0], pca_features2_sol0[1], label = 'low solubility')
plt.scatter(pca_features2_sol1[0], pca_features2_sol1[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2[0], pca_features2_sol2[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for training and validation dataset split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_split_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features1_samples[0], pca_features1_samples[1], label = 'no condition')
plt.scatter(pca_features1_sol0_samples[0], pca_features1_sol0_samples[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_samples[0], pca_features1_sol1_samples[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_samples[0], pca_features1_sol2_samples[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated samples split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_samples_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features2_samples[0], pca_features2_samples[1], label = 'no condition')
plt.scatter(pca_features2_sol0_samples[0], pca_features2_sol0_samples[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_samples[0], pca_features2_sol1_samples[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_samples[0], pca_features2_sol2_samples[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated samples split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_samples_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features1_variants[0], pca_features1_variants[1], label = 'no condition')
plt.scatter(pca_features1_sol0_variants[0], pca_features1_sol0_variants[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_variants[0], pca_features1_sol1_variants[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_variants[0], pca_features1_sol2_variants[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_variants_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features2_variants[0], pca_features2_variants[1], label = 'no condition')
plt.scatter(pca_features2_sol0_variants[0], pca_features2_sol0_variants[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_variants[0], pca_features2_sol1_variants[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_variants[0], pca_features2_sol2_variants[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_variants_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features1_cond[0], pca_features1_cond[1], label = 'train and validation data')
plt.scatter(pca_features1_sol0_cond[0], pca_features1_sol0_cond[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_cond[0], pca_features1_sol1_cond[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_cond[0], pca_features1_sol2_cond[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated conditional variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_cond_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2_cond[0], pca_features2_cond[1], label = 'train and validation data')
plt.scatter(pca_features2_sol0_cond[0], pca_features2_sol0_cond[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_cond[0], pca_features2_sol1_cond[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_cond[0], pca_features2_sol2_cond[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated conditional variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/msa_cond_PCA2.png', bbox_inches='tight') 
plt.close()

outputs_aaletters = np.load('../../output/generated_PCA/ar_outputs_aaletters.npy')

outputs_aaletters_sol0 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol0.npy')
outputs_aaletters_sol1 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol1.npy')
outputs_aaletters_sol2 = np.load('../../output/generated_PCA/ar_outputs_aaletters_sol2.npy') 

outputs_aaletters_cond_sol0 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol0.npy')
outputs_aaletters_cond_sol1 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol1.npy')
outputs_aaletters_cond_sol2 = np.load('../../output/generated_PCA/ar_outputs_aaletters_cond_sol2.npy') 

outputs = np.load('../../output/generated_PCA/ar_outputs.npy')
outputs_samples = np.load('../../output/generated_PCA/ar_outputs_samples.npy')
outputs_variants = np.load('../../output/generated_PCA/ar_outputs_variants.npy')

outputs_sol0 = np.load('../../output/generated_PCA/ar_outputs_sol0.npy')
outputs_sol0_samples = np.load('../../output/generated_PCA/ar_outputs_sol0_samples.npy')
outputs_sol0_variants = np.load('../../output/generated_PCA/ar_outputs_sol0_variants.npy')

outputs_sol1 = np.load('../../output/generated_PCA/ar_outputs_sol1.npy')
outputs_sol1_samples = np.load('../../output/generated_PCA/ar_outputs_sol1_samples.npy')
outputs_sol1_variants = np.load('../../output/generated_PCA/ar_outputs_sol1_variants.npy')

outputs_sol2 = np.load('../../output/generated_PCA/ar_outputs_sol2.npy')
outputs_sol2_samples = np.load('../../output/generated_PCA/ar_outputs_sol2_samples.npy')
outputs_sol2_variants = np.load('../../output/generated_PCA/ar_outputs_sol2_variants.npy')

outputs_cond = np.load('../../output/generated_PCA/ar_outputs_cond.npy')
outputs_cond_sol0 = np.load('../../output/generated_PCA/ar_outputs_cond_sol0.npy')
outputs_cond_sol1 = np.load('../../output/generated_PCA/ar_outputs_cond_sol1.npy')
outputs_cond_sol2 = np.load('../../output/generated_PCA/ar_outputs_cond_sol2.npy')

pca_features1_aaletters = pca.fit_transform(outputs_aaletters[0]).reshape(2, -1)
pca_features2_aaletters = pca.fit_transform(outputs_aaletters[1]).reshape(2, -1)

pca_features1_aaletters_sol0 = pca.fit_transform(outputs_aaletters_sol0[0]).reshape(2, -1)
pca_features2_aaletters_sol0 = pca.fit_transform(outputs_aaletters_sol0[1]).reshape(2, -1)

pca_features1_aaletters_sol1 = pca.fit_transform(outputs_aaletters_sol1[0]).reshape(2, -1)
pca_features2_aaletters_sol1 = pca.fit_transform(outputs_aaletters_sol1[1]).reshape(2, -1)

pca_features1_aaletters_sol2 = pca.fit_transform(outputs_aaletters_sol2[0]).reshape(2, -1)
pca_features2_aaletters_sol2 = pca.fit_transform(outputs_aaletters_sol2[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol0 = pca.fit_transform(outputs_aaletters_cond_sol0[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol0 = pca.fit_transform(outputs_aaletters_cond_sol0[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol1 = pca.fit_transform(outputs_aaletters_cond_sol1[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol1 = pca.fit_transform(outputs_aaletters_cond_sol1[1]).reshape(2, -1)

pca_features1_aaletters_cond_sol2 = pca.fit_transform(outputs_aaletters_cond_sol2[0]).reshape(2, -1)
pca_features2_aaletters_cond_sol2 = pca.fit_transform(outputs_aaletters_cond_sol2[1]).reshape(2, -1)

pca_features1 = pca.fit_transform(outputs[0]).reshape(2, -1)
pca_features2 = pca.fit_transform(outputs[1]).reshape(2, -1)

pca_features1_samples = pca.fit_transform(outputs_samples[0]).reshape(2, -1)
pca_features2_samples = pca.fit_transform(outputs_samples[1]).reshape(2, -1)

pca_features1_variants = pca.fit_transform(outputs_variants[0]).reshape(2, -1)
pca_features2_variants = pca.fit_transform(outputs_variants[1]).reshape(2, -1)

pca_features1_sol0 = pca.fit_transform(outputs_sol0[0]).reshape(2, -1)
pca_features2_sol0 = pca.fit_transform(outputs_sol0[1]).reshape(2, -1)

pca_features1_sol0_samples = pca.fit_transform(outputs_sol0_samples[0]).reshape(2, -1)
pca_features2_sol0_samples = pca.fit_transform(outputs_sol0_samples[1]).reshape(2, -1)

pca_features1_sol0_variants = pca.fit_transform(outputs_sol0_variants[0]).reshape(2, -1)
pca_features2_sol0_variants = pca.fit_transform(outputs_sol0_variants[1]).reshape(2, -1)

pca_features1_sol1 = pca.fit_transform(outputs_sol1[0]).reshape(2, -1)
pca_features2_sol1 = pca.fit_transform(outputs_sol1[1]).reshape(2, -1)

pca_features1_sol1_samples = pca.fit_transform(outputs_sol1_samples[0]).reshape(2, -1)
pca_features2_sol1_samples = pca.fit_transform(outputs_sol1_samples[1]).reshape(2, -1)

pca_features1_sol1_variants = pca.fit_transform(outputs_sol1_variants[0]).reshape(2, -1)
pca_features2_sol1_variants = pca.fit_transform(outputs_sol1_variants[1]).reshape(2, -1)

pca_features1_sol2 = pca.fit_transform(outputs_sol2[0]).reshape(2, -1)
pca_features2_sol2 = pca.fit_transform(outputs_sol2[1]).reshape(2, -1)

pca_features1_sol2_samples = pca.fit_transform(outputs_sol2_samples[0]).reshape(2, -1)
pca_features2_sol2_samples = pca.fit_transform(outputs_sol2_samples[1]).reshape(2, -1)

pca_features1_sol2_variants = pca.fit_transform(outputs_sol2_variants[0]).reshape(2, -1)
pca_features2_sol2_variants = pca.fit_transform(outputs_sol2_variants[1]).reshape(2, -1)

pca_features1_cond = pca.fit_transform(outputs_cond[0]).reshape(2, -1)
pca_features2_cond = pca.fit_transform(outputs_cond[1]).reshape(2, -1)

pca_features1_sol0_cond = pca.fit_transform(outputs_cond_sol0[0]).reshape(2, -1)
pca_features2_sol0_cond = pca.fit_transform(outputs_cond_sol0[1]).reshape(2, -1)

pca_features1_sol1_cond = pca.fit_transform(outputs_cond_sol1[0]).reshape(2, -1)
pca_features2_sol1_cond = pca.fit_transform(outputs_cond_sol1[1]).reshape(2, -1)

pca_features1_sol2_cond = pca.fit_transform(outputs_cond_sol2[0]).reshape(2, -1)
pca_features2_sol2_cond = pca.fit_transform(outputs_cond_sol2[1]).reshape(2, -1)

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters[0][letter_index], pca_features1_aaletters[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters[0][letter_index], pca_features1_aaletters[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters[0][letter_index], pca_features2_aaletters[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters[0][letter_index], pca_features2_aaletters[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol0[0][letter_index], pca_features1_aaletters_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol0[0][letter_index], pca_features1_aaletters_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol0_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol0[0][letter_index], pca_features2_aaletters_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol0[0][letter_index], pca_features2_aaletters_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol0_PCA2.png', bbox_inches='tight') 
plt.close()
 
for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol1[0][letter_index], pca_features1_aaletters_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol1[0][letter_index], pca_features1_aaletters_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol1_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol1[0][letter_index], pca_features2_aaletters_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol1[0][letter_index], pca_features2_aaletters_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol1_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_sol2[0][letter_index], pca_features1_aaletters_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_sol2[0][letter_index], pca_features1_aaletters_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol2_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_sol2[0][letter_index], pca_features2_aaletters_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_sol2[0][letter_index], pca_features2_aaletters_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_sol2_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol0[0][letter_index], pca_features1_aaletters_cond_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol0[0][letter_index], pca_features1_aaletters_cond_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol0_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol0[0][letter_index], pca_features2_aaletters_cond_sol0[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol0[0][letter_index], pca_features2_aaletters_cond_sol0[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with low solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol0_PCA2.png', bbox_inches='tight') 
plt.close()
 
for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol1[0][letter_index], pca_features1_aaletters_cond_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol1[0][letter_index], pca_features1_aaletters_cond_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol1_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol1[0][letter_index], pca_features2_aaletters_cond_sol1[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol1[0][letter_index], pca_features2_aaletters_cond_sol1[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with mid solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol1_PCA2.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features1_aaletters_cond_sol2[0][letter_index], pca_features1_aaletters_cond_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features1_aaletters_cond_sol2[0][letter_index], pca_features1_aaletters_cond_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol2_PCA1.png', bbox_inches='tight') 
plt.close()

for letter_index in range(1, len(aaletters)):
	plot_letter(pca_features2_aaletters_cond_sol2[0][letter_index], pca_features2_aaletters_cond_sol2[1][letter_index], aaletters[letter_index])
	plt.annotate(aaletters[letter_index], (pca_features2_aaletters_cond_sol2[0][letter_index], pca_features2_aaletters_cond_sol2[1][letter_index]))
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for aminoacids with high solubility model")
plt.legend(ncol = 3, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_aaletters_cond_sol2_PCA2.png', bbox_inches='tight') 
plt.close()
 
plt.scatter(pca_features1[0], pca_features1[1], label = 'all data')
plt.scatter(pca_features1_sol0[0], pca_features1_sol0[1], label = 'low solubility')
plt.scatter(pca_features1_sol1[0], pca_features1_sol1[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2[0], pca_features1_sol2[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for training and validation dataset split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_split_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'all data')
plt.scatter(pca_features2_sol0[0], pca_features2_sol0[1], label = 'low solubility')
plt.scatter(pca_features2_sol1[0], pca_features2_sol1[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2[0], pca_features2_sol2[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for training and validation dataset split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_split_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features1_samples[0], pca_features1_samples[1], label = 'no condition')
plt.scatter(pca_features1_sol0_samples[0], pca_features1_sol0_samples[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_samples[0], pca_features1_sol1_samples[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_samples[0], pca_features1_sol2_samples[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated samples split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_samples_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features2_samples[0], pca_features2_samples[1], label = 'no condition')
plt.scatter(pca_features2_sol0_samples[0], pca_features2_sol0_samples[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_samples[0], pca_features2_sol1_samples[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_samples[0], pca_features2_sol2_samples[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated samples split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_samples_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features1_variants[0], pca_features1_variants[1], label = 'no condition')
plt.scatter(pca_features1_sol0_variants[0], pca_features1_sol0_variants[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_variants[0], pca_features1_sol1_variants[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_variants[0], pca_features1_sol2_variants[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_variants_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2[0], pca_features2[1], label = 'train and validation data')
plt.scatter(pca_features2_variants[0], pca_features2_variants[1], label = 'no condition')
plt.scatter(pca_features2_sol0_variants[0], pca_features2_sol0_variants[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_variants[0], pca_features2_sol1_variants[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_variants[0], pca_features2_sol2_variants[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_variants_PCA2.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features1_cond[0], pca_features1_cond[1], label = 'train and validation data')
plt.scatter(pca_features1_sol0_cond[0], pca_features1_sol0_cond[1], label = 'low solubility')
plt.scatter(pca_features1_sol1_cond[0], pca_features1_sol1_cond[1], label = 'mid solubility')
plt.scatter(pca_features1_sol2_cond[0], pca_features1_sol2_cond[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated conditional variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_cond_PCA1.png', bbox_inches='tight') 
plt.close()

plt.scatter(pca_features2_cond[0], pca_features2_cond[1], label = 'train and validation data')
plt.scatter(pca_features2_sol0_cond[0], pca_features2_sol0_cond[1], label = 'low solubility')
plt.scatter(pca_features2_sol1_cond[0], pca_features2_sol1_cond[1], label = 'mid solubility')
plt.scatter(pca_features2_sol2_cond[0], pca_features2_sol2_cond[1], label = 'high solubility')
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title("PCA for generated conditional variants split by solubility")
plt.legend(ncol = 2, bbox_to_anchor=(0, -0.4), loc="lower left")
plt.savefig('../results/read_PCA/ar_cond_PCA2.png', bbox_inches='tight') 
plt.close()