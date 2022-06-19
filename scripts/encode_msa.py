import numpy as np
from models.vaes import MSAVAE
from utils.io import load_gzdata, read_fasta
from utils.data_loaders import to_one_hot, right_pad 
from utils.conditions_loading import load_conditions
from utils import aa_letters
aaletters = aa_letters

# Load aligned sequences

names, msa_seqs = load_gzdata('data/training_data/luxafilt_llmsa_train.fa.gz', one_hot=False)
_, val_msa_seqs = load_gzdata('data/training_data/luxafilt_llmsa_val.fa.gz', one_hot=False)

conditions_train = load_conditions('data/training_data/conditions_ll_train.txt') 
conditions_val = load_conditions('data/training_data/conditions_ll_val.txt') 

_, msa_seqs_sol0 = load_gzdata('data/training_data/luxafilt_llmsa_train_sol0.fa.gz', one_hot=False)
_, val_msa_seqs_sol0 = load_gzdata('data/training_data/luxafilt_llmsa_val_sol0.fa.gz', one_hot=False)

_, msa_seqs_sol1 = load_gzdata('data/training_data/luxafilt_llmsa_train_sol1.fa.gz', one_hot=False)
_, val_msa_seqs_sol1 = load_gzdata('data/training_data/luxafilt_llmsa_val_sol1.fa.gz', one_hot=False)

_, msa_seqs_sol2 = load_gzdata('data/training_data/luxafilt_llmsa_train_sol2.fa.gz', one_hot=False)
_, val_msa_seqs_sol2 = load_gzdata('data/training_data/luxafilt_llmsa_val_sol2.fa.gz', one_hot=False)

_, msa_variants = read_fasta('../output/generated_sequences/msavae_variants.fa')
_, msa_variants_sol0 = read_fasta('../output/generated_sequences/msavae_sol0_variants.fa')
_, msa_variants_sol1 = read_fasta('../output/generated_sequences/msavae_sol1_variants.fa')
_, msa_variants_sol2 = read_fasta('../output/generated_sequences/msavae_sol2_variants.fa')

_, msa_samples = read_fasta('../output/generated_sequences/msavae_samples.fa')
_, msa_samples_sol0 = read_fasta('../output/generated_sequences/msavae_sol0_samples.fa')
_, msa_samples_sol1 = read_fasta('../output/generated_sequences/msavae_sol1_samples.fa')
_, msa_samples_sol2 = read_fasta('../output/generated_sequences/msavae_sol2_samples.fa')

_, msa_with_conditions_variants_sol0 = read_fasta('../output/generated_sequences/msavae_with_conditions_sol0_variants.fa')
_, msa_with_conditions_variants_sol1 = read_fasta('../output/generated_sequences/msavae_with_conditions_sol1_variants.fa')
_, msa_with_conditions_variants_sol2 = read_fasta('../output/generated_sequences/msavae_with_conditions_sol2_variants.fa')

target_len_for_seq = 360

# Define data

train_gen_aaletters = to_one_hot(right_pad(aaletters, target_len_for_seq))
train_gen = to_one_hot(right_pad(msa_seqs + val_msa_seqs, target_len_for_seq))

train_gen_sol0 = to_one_hot(right_pad(msa_seqs_sol0 + val_msa_seqs_sol0, target_len_for_seq))
train_gen_sol1 = to_one_hot(right_pad(msa_seqs_sol1 + val_msa_seqs_sol1, target_len_for_seq))
train_gen_sol2 = to_one_hot(right_pad(msa_seqs_sol2 + val_msa_seqs_sol2, target_len_for_seq))

train_gen_samples = to_one_hot(right_pad(msa_samples, target_len_for_seq))
train_gen_sol0_samples = to_one_hot(right_pad(msa_samples_sol0, target_len_for_seq))
train_gen_sol1_samples = to_one_hot(right_pad(msa_samples_sol1, target_len_for_seq))
train_gen_sol2_samples = to_one_hot(right_pad(msa_samples_sol2, target_len_for_seq))

train_gen_variants = to_one_hot(right_pad(msa_variants, target_len_for_seq))
train_gen_sol0_variants = to_one_hot(right_pad(msa_variants_sol0, target_len_for_seq))
train_gen_sol1_variants = to_one_hot(right_pad(msa_variants_sol1, target_len_for_seq))
train_gen_sol2_variants = to_one_hot(right_pad(msa_variants_sol2, target_len_for_seq))

# Define model
print('Building model')

model = MSAVAE(original_dim=360, latent_dim=10) 
model_sol0 = MSAVAE(original_dim=360, latent_dim=10)
model_sol1 = MSAVAE(original_dim=360, latent_dim=10)
model_sol2 = MSAVAE(original_dim=360, latent_dim=10)
model_cond = MSAVAE(original_dim=360, latent_dim=10, n_conditions=3)
 
model.load_weights('data/weights/msavae.h5') 
model_sol0.load_weights('../output/weights/msavae_sol0.h5')
model_sol1.load_weights('../output/weights/msavae_sol1.h5')
model_sol2.load_weights('../output/weights/msavae_sol2.h5')
model_cond.load_weights('../output/weights/msavae_with_conditions.h5')

outputs_cond = model_cond.E.predict([train_gen, np.array(conditions_train + conditions_val)])
np.save('../output/generated_PCA/msa_outputs_cond.npy', outputs_cond)
 
sol0 = []
for i in range(len(msa_with_conditions_variants_sol0)):
    sol0.append(np.array([1,0,0]))
sol0 = np.array(sol0)
outputs_cond_sol0 = model_cond.E.predict([train_gen_sol0_variants, sol0])
np.save('../output/generated_PCA/msa_outputs_cond_sol0.npy', outputs_cond_sol0)
sol0 = []
for i in range(len(aaletters)):
    sol0.append(np.array([1,0,0]))
sol0 = np.array(sol0)
outputs_aaletters_cond_sol0 = model_cond.E.predict([train_gen_aaletters, sol0])
np.save('../output/generated_PCA/msa_outputs_aaletters_cond_sol0.npy', outputs_aaletters_cond_sol0)

sol1 = []
for i in range(len(msa_with_conditions_variants_sol1)):
    sol1.append(np.array([0,1,0]))
sol1 = np.array(sol1)
outputs_cond_sol1 = model_cond.E.predict([train_gen_sol1_variants, sol1])
np.save('../output/generated_PCA/msa_outputs_cond_sol1.npy', outputs_cond_sol1)
sol1 = []
for i in range(len(aaletters)):
    sol1.append(np.array([0,1,0]))
sol1 = np.array(sol1)
outputs_aaletters_cond_sol1 = model_cond.E.predict([train_gen_aaletters, sol1])
np.save('../output/generated_PCA/msa_outputs_aaletters_cond_sol1.npy', outputs_aaletters_cond_sol1)

sol2 = []
for i in range(len(msa_with_conditions_variants_sol2)):
    sol2.append(np.array([0,0,1]))
sol2 = np.array(sol2)
outputs_cond_sol2 = model_cond.E.predict([train_gen_sol2_variants, sol2])
np.save('../output/generated_PCA/msa_outputs_cond_sol2.npy', outputs_cond_sol2)
sol2 = []
for i in range(len(aaletters)):
    sol2.append(np.array([0,0,1]))
sol2 = np.array(sol2)
outputs_aaletters_cond_sol2 = model_cond.E.predict([train_gen_aaletters, sol2])
np.save('../output/generated_PCA/msa_outputs_aaletters_cond_sol2.npy', outputs_aaletters_cond_sol2)

outputs_aaletters = model.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/msa_outputs_aaletters.npy', outputs_aaletters)
outputs_aaletters_sol0 = model_sol0.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/msa_outputs_aaletters_sol0.npy', outputs_aaletters_sol0)
outputs_aaletters_sol1 = model_sol1.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/msa_outputs_aaletters_sol1.npy', outputs_aaletters_sol1)
outputs_aaletters_sol2 = model_sol2.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/msa_outputs_aaletters_sol2.npy', outputs_aaletters_sol2)

outputs = model.E.predict(train_gen)
np.save('../output/generated_PCA/msa_outputs.npy', outputs)
outputs_samples = model.E.predict(train_gen_samples)
np.save('../output/generated_PCA/msa_outputs_samples.npy', outputs_samples)
outputs_variants = model.E.predict(train_gen_variants)
np.save('../output/generated_PCA/msa_outputs_variants.npy', outputs_variants)

outputs_sol0 = model_sol0.E.predict(train_gen_sol0)
np.save('../output/generated_PCA/msa_outputs_sol0.npy', outputs_sol0)
outputs_sol0_samples = model_sol0.E.predict(train_gen_sol0_samples)
np.save('../output/generated_PCA/msa_outputs_sol0_samples.npy', outputs_sol0_samples)
outputs_sol0_variants = model_sol0.E.predict(train_gen_sol0_variants)
np.save('../output/generated_PCA/msa_outputs_sol0_variants.npy', outputs_sol0_variants)

outputs_sol1 = model_sol1.E.predict(train_gen_sol1)
np.save('../output/generated_PCA/msa_outputs_sol1.npy', outputs_sol1)
outputs_sol1_samples = model_sol1.E.predict(train_gen_sol1_samples)
np.save('../output/generated_PCA/msa_outputs_sol1_samples.npy', outputs_sol1_samples)
outputs_sol1_variants = model_sol1.E.predict(train_gen_sol1_variants)
np.save('../output/generated_PCA/msa_outputs_sol1_variants.npy', outputs_sol1_variants)

outputs_sol2 = model_sol2.E.predict(train_gen_sol2)
np.save('../output/generated_PCA/msa_outputs_sol2.npy', outputs_sol2)
outputs_sol2_samples = model_sol2.E.predict(train_gen_sol2_samples)
np.save('../output/generated_PCA/msa_outputs_sol2_samples.npy', outputs_sol2_samples)
outputs_sol2_variants = model_sol2.E.predict(train_gen_sol2_variants)
np.save('../output/generated_PCA/msa_outputs_sol2_variants.npy', outputs_sol2_variants)

print("Finished")