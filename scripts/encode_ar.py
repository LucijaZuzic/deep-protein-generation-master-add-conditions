import numpy as np
from models.vaes import ARVAE
from utils.io import load_gzdata, read_fasta
from utils.data_loaders import to_one_hot, right_pad 
from utils.conditions_loading import load_conditions
from utils import aa_letters
aaletters = aa_letters

# Load aligned sequences

names, ar_seqs = load_gzdata('data/training_data/ll_train.fa.gz', one_hot=False)
_, val_ar_seqs = load_gzdata('data/training_data/ll_val.fa.gz', one_hot=False)

conditions_train = load_conditions('data/training_data/conditions_ll_train.txt') 
conditions_val = load_conditions('data/training_data/conditions_ll_val.txt') 

_, ar_seqs_sol0 = load_gzdata('data/training_data/ll_train_sol0.fa.gz', one_hot=False)
_, val_ar_seqs_sol0 = load_gzdata('data/training_data/ll_val_sol0.fa.gz', one_hot=False)

_, ar_seqs_sol1 = load_gzdata('data/training_data/ll_train_sol1.fa.gz', one_hot=False)
_, val_ar_seqs_sol1 = load_gzdata('data/training_data/ll_val_sol1.fa.gz', one_hot=False)

_, ar_seqs_sol2 = load_gzdata('data/training_data/ll_train_sol2.fa.gz', one_hot=False)
_, val_ar_seqs_sol2 = load_gzdata('data/training_data/ll_val_sol2.fa.gz', one_hot=False)

_, ar_variants = read_fasta('../output/generated_sequences/arvae_variants.fa')
_, ar_variants_sol0 = read_fasta('../output/generated_sequences/arvae_sol0_variants.fa')
_, ar_variants_sol1 = read_fasta('../output/generated_sequences/arvae_sol1_variants.fa')
_, ar_variants_sol2 = read_fasta('../output/generated_sequences/arvae_sol2_variants.fa')

_, ar_samples = read_fasta('../output/generated_sequences/arvae_samples.fa')
_, ar_samples_sol0 = read_fasta('../output/generated_sequences/arvae_sol0_samples.fa')
_, ar_samples_sol1 = read_fasta('../output/generated_sequences/arvae_sol1_samples.fa')
_, ar_samples_sol2 = read_fasta('../output/generated_sequences/arvae_sol2_samples.fa')

_, ar_with_conditions_variants_sol0 = read_fasta('../output/generated_sequences/arvae_with_conditions_sol0_variants.fa')
_, ar_with_conditions_variants_sol1 = read_fasta('../output/generated_sequences/arvae_with_conditions_sol1_variants.fa')
_, ar_with_conditions_variants_sol2 = read_fasta('../output/generated_sequences/arvae_with_conditions_sol2_variants.fa')

target_len_for_seq = 504

# Define data

train_gen_aaletters = to_one_hot(right_pad(aaletters, target_len_for_seq))
train_gen = to_one_hot(right_pad(ar_seqs + val_ar_seqs, target_len_for_seq))

train_gen_sol0 = to_one_hot(right_pad(ar_seqs_sol0 + val_ar_seqs_sol0, target_len_for_seq))
train_gen_sol1 = to_one_hot(right_pad(ar_seqs_sol1 + val_ar_seqs_sol1, target_len_for_seq))
train_gen_sol2 = to_one_hot(right_pad(ar_seqs_sol2 + val_ar_seqs_sol2, target_len_for_seq))

train_gen_samples = to_one_hot(right_pad(ar_samples, target_len_for_seq))
train_gen_sol0_samples = to_one_hot(right_pad(ar_samples_sol0, target_len_for_seq))
train_gen_sol1_samples = to_one_hot(right_pad(ar_samples_sol1, target_len_for_seq))
train_gen_sol2_samples = to_one_hot(right_pad(ar_samples_sol2, target_len_for_seq))

train_gen_variants = to_one_hot(right_pad(ar_variants, target_len_for_seq))
train_gen_sol0_variants = to_one_hot(right_pad(ar_variants_sol0, target_len_for_seq))
train_gen_sol1_variants = to_one_hot(right_pad(ar_variants_sol1, target_len_for_seq))
train_gen_sol2_variants = to_one_hot(right_pad(ar_variants_sol2, target_len_for_seq))

# Define model
print('Building model')

model = ARVAE() 
model_sol0 = ARVAE()
model_sol1 = ARVAE()
model_sol2 = ARVAE()
model_cond = ARVAE(n_conditions=3)
 
model.load_weights('data/weights/arvae.h5') 
model_sol0.load_weights('../output/weights/arvae_sol0.h5')
model_sol1.load_weights('../output/weights/arvae_sol1.h5')
model_sol2.load_weights('../output/weights/arvae_sol2.h5')
model_cond.load_weights('../output/weights/arvae_with_conditions.h5')

print('Building model cond')
outputs_cond = model_cond.E.predict([train_gen, np.array(conditions_train + conditions_val)])
np.save('../output/generated_PCA/ar_outputs_cond.npy', outputs_cond)
 
print('Building model cond sol0')
sol0 = []
for i in range(len(ar_with_conditions_variants_sol0)):
    sol0.append(np.array([1,0,0]))
sol0 = np.array(sol0)
outputs_cond_sol0 = model_cond.E.predict([train_gen_sol0_variants, sol0])
np.save('../output/generated_PCA/ar_outputs_cond_sol0.npy', outputs_cond_sol0)
sol0 = []
for i in range(len(aaletters)):
    sol0.append(np.array([1,0,0]))
sol0 = np.array(sol0)
outputs_aaletters_cond_sol0 = model_cond.E.predict([train_gen_aaletters, sol0])
np.save('../output/generated_PCA/ar_outputs_aaletters_cond_sol0.npy', outputs_aaletters_cond_sol0)

print('Building model cond sol1')
sol1 = []
for i in range(len(ar_with_conditions_variants_sol1)):
    sol1.append(np.array([0,1,0]))
sol1 = np.array(sol1)
outputs_cond_sol1 = model_cond.E.predict([train_gen_sol1_variants, sol1])
np.save('../output/generated_PCA/ar_outputs_cond_sol1.npy', outputs_cond_sol1)
sol1 = []
for i in range(len(aaletters)):
    sol1.append(np.array([0,1,0]))
sol1 = np.array(sol1)
outputs_aaletters_cond_sol1 = model_cond.E.predict([train_gen_aaletters, sol1])
np.save('../output/generated_PCA/ar_outputs_aaletters_cond_sol1.npy', outputs_aaletters_cond_sol1)

print('Building model cond sol2')
sol2 = []
for i in range(len(ar_with_conditions_variants_sol2)):
    sol2.append(np.array([0,0,1]))
sol2 = np.array(sol2)
outputs_cond_sol2 = model_cond.E.predict([train_gen_sol2_variants, sol2])
np.save('../output/generated_PCA/ar_outputs_cond_sol2.npy', outputs_cond_sol2)
sol2 = []
for i in range(len(aaletters)):
    sol2.append(np.array([0,0,1]))
sol2 = np.array(sol2)
outputs_aaletters_cond_sol2 = model_cond.E.predict([train_gen_aaletters, sol2])
np.save('../output/generated_PCA/ar_outputs_aaletters_cond_sol2.npy', outputs_aaletters_cond_sol2)

print('Building model letters')
outputs_aaletters = model.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/ar_outputs_aaletters.npy', outputs_aaletters)
print('Building model letters sol0')
outputs_aaletters_sol0 = model_sol0.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/ar_outputs_aaletters_sol0.npy', outputs_aaletters_sol0)
print('Building model letters sol1')
outputs_aaletters_sol1 = model_sol1.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/ar_outputs_aaletters_sol1.npy', outputs_aaletters_sol1)
print('Building model letters sol2')
outputs_aaletters_sol2 = model_sol2.E.predict(train_gen_aaletters)
np.save('../output/generated_PCA/ar_outputs_aaletters_sol2.npy', outputs_aaletters_sol2)

print('Building model all')
outputs = model.E.predict(train_gen)
np.save('../output/generated_PCA/ar_outputs.npy', outputs)
print('Building model samples')
outputs_samples = model.E.predict(train_gen_samples)
np.save('../output/generated_PCA/ar_outputs_samples.npy', outputs_samples)
print('Building model variants')
outputs_variants = model.E.predict(train_gen_variants)
np.save('../output/generated_PCA/ar_outputs_variants.npy', outputs_variants)

print('Building model all sol0')
outputs_sol0 = model_sol0.E.predict(train_gen_sol0)
np.save('../output/generated_PCA/ar_outputs_sol0.npy', outputs_sol0)
print('Building model samples sol0')
outputs_sol0_samples = model_sol0.E.predict(train_gen_sol0_samples)
np.save('../output/generated_PCA/ar_outputs_sol0_samples.npy', outputs_sol0_samples)
print('Building model variants sol0')
outputs_sol0_variants = model_sol0.E.predict(train_gen_sol0_variants)
np.save('../output/generated_PCA/ar_outputs_sol0_variants.npy', outputs_sol0_variants)

print('Building model all sol1')
outputs_sol1 = model_sol1.E.predict(train_gen_sol1)
np.save('../output/generated_PCA/ar_outputs_sol1.npy', outputs_sol1)
print('Building model samples sol1')
outputs_sol1_samples = model_sol1.E.predict(train_gen_sol1_samples)
np.save('../output/generated_PCA/ar_outputs_sol1_samples.npy', outputs_sol1_samples)
print('Building model variants sol1')
outputs_sol1_variants = model_sol1.E.predict(train_gen_sol1_variants)
np.save('../output/generated_PCA/ar_outputs_sol1_variants.npy', outputs_sol1_variants)

print('Building model all sol2')
outputs_sol2 = model_sol2.E.predict(train_gen_sol2)
np.save('../output/generated_PCA/ar_outputs_sol2.npy', outputs_sol2)
print('Building model samples sol2')
outputs_sol2_samples = model_sol2.E.predict(train_gen_sol2_samples)
np.save('../output/generated_PCA/ar_outputs_sol2_samples.npy', outputs_sol2_samples)
print('Building model variants sol2')
outputs_sol2_variants = model_sol2.E.predict(train_gen_sol2_variants)
np.save('../output/generated_PCA/ar_outputs_sol2_variants.npy', outputs_sol2_variants)

print("Finished")