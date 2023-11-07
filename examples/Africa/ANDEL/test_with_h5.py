# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 21:02:02 2023

@author: dboateng
"""

import numpy as np
import h5py
from sklearn.model_selection import train_test_split

def train_test(file):
    
    with h5py.File(file, 'r') as h5:
        X_train = h5['X_train'][:,:,256:,:,:]
        X_valid = h5['X_valid'][:,:,256:,:,:]
        y_train = h5['Y_train'][:,:,256:,:,:]
        y_valid = h5['Y_valid'][:,:,256:,:,:]

    return X_train, X_valid, y_train, y_valid


path_to_file = "D:/Datasets/IMERG/practice_data/complete_final_dataset_IMERG_Nov_05.h5"

X_t, X_v, y_t, y_v = train_test(path_to_file)
