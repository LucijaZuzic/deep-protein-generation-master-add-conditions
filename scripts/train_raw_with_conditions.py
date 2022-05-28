import tensorflow as tf

import argparse
import numpy as np
from keras.callbacks import CSVLogger, ModelCheckpoint

from models.vaes import ARVAE
from utils.io import load_gzdata
from utils.data_loaders import one_hot_generator
from utils.conditions_loading import load_conditions

# Define training parameters
batch_size = 32
seed = 0
n_epochs = 50
verbose = 1
save_all_epochs = False

seed = np.random.seed(seed)

# Load aligned sequences
_, raw_seqs = load_gzdata('data/training_data/ll_train.fa.gz', one_hot=False)
_, val_raw_seqs = load_gzdata('data/training_data/ll_val.fa.gz', one_hot=False)

# Define data generators

conditions_train = load_conditions('data/training_data/conditions_ll_train.txt') 
conditions_val = load_conditions('data/training_data/conditions_ll_val.txt') 

train_gen = one_hot_generator(raw_seqs, padding=504, conditions = conditions_train)
val_gen = one_hot_generator(val_raw_seqs, padding=504, conditions = conditions_val)

# Define model
print('Building model')
model = ARVAE(n_conditions=3)

# (Optionally) define callbacks
callbacks=[CSVLogger('output/logs/arvae_with_conditions.csv')]

if save_all_epochs:
    callbacks.append(ModelCheckpoint('output/weights/arvae_with_conditions'+'.{epoch:02d}-{luxa_errors_mean:.2f}.hdf5',
                                     save_best_only=False, verbose=1))

# Train model https://github.com/keras-team/keras/issues/8595
model.VAE.fit_generator(generator=train_gen,
                        steps_per_epoch=len(raw_seqs) // batch_size,
                        verbose=verbose,
                        epochs=n_epochs,
                        validation_data=val_gen,
                        validation_steps=len(val_raw_seqs) // batch_size,
                        callbacks=callbacks)

if not save_all_epochs:
  model.save_weights('output/weights/arvae_with_conditions.h5')