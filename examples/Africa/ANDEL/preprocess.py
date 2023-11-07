# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 18:32:43 2023

@author: dboateng
"""

import os
import numpy as np
import pandas as pd 
import xarray as xr
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

import pickle


path_to_data = "D:/Datasets/IMERG/"

datasets = xr.open_mfdataset(path_to_data+"*.nc4")
def train_test():
    prec_data = datasets.precipitationCal.compute()
    data = np.expand_dims(prec_data, axis=-1)
    
    indx = np.arange(0, data.shape[0], 9+3)
    X = [data[indx[i]:indx[i+1] - 3, :, :, 0] for i in range(len(indx)-1)]
    X = np.expand_dims(X, axis=-1)
    
    Y = [data[indx[i]+ 9:indx[i+1], :, :, 0] - data[indx[i+1]- 3-1, :, :, 0]  for i in range(len(indx)-1)]
    Y =  np.expand_dims(Y, axis=-1)
    Y =  np.expand_dims(Y, axis=-1)
    X_train, X_valid, y_train, y_valid = train_test_split(X, Y, test_size=0.2, random_state=42)
    
    
    percentage = np.array([(X_train[i, 0].squeeze()>1).mean()*100 for i in range(X_train.shape[0])])
    percentage_valid = np.array([(X_valid[i, 0].squeeze()>1).mean()*100 for i in range(X_valid.shape[0])])
    
    print(X_train.shape, y_valid.shape)
    
    return X_train[percentage>1], X_valid[percentage_valid>1],y_train[percentage>1], y_valid[percentage_valid>1]


X_train, X_valid, y_train, y_valid = train_test()

def reshape(arr):
    arr = arr.squeeze()
    arr = arr.reshape(arr.shape[0],arr.shape[1],arr.shape[2]*arr.shape[3])
    arr = arr.swapaxes(0,1)
    arr = arr.reshape(arr.shape[0],arr.shape[1]*arr.shape[2]).T

    return arr


X_train = reshape(X_train)
y_train = reshape(y_train)
# X_valid = reshape(X_train)
# y_valid = reshape(y_train)

class LinearRegression_Class:

    def load_model(path):
        mdl = LinearRegression_Class()
        with open(path,'rb') as f:
            mdl.sk_model = pickle.load(f)

        return mdl

    def reshape(arr):
        arr = arr.squeeze()
        arr = arr.reshape(arr.shape[0],arr.shape[1]*arr.shape[2])
        arr = arr.swapaxes(0,1)
        return arr

    def reshape_back(arr,s):
        arr = arr.swapaxes(0,1)
        arr = arr.reshape(s[0],s[1],s[2])
        return arr

    def fit(self, X_train, y_train):
        self.sk_model = LinearRegression().fit(X_train, y_train)

    def predict(self,X):
        s = X.squeeze().shape
        
        reshaped_X = LinearRegression_Class.reshape(X)
        reshaped_Y = self.sk_model.predict(reshaped_X)
        Y = LinearRegression_Class.reshape_back(reshaped_Y,(3,s[1],s[2]))

        return Y
    
mdl = LinearRegression_Class()
mdl.fit(X_train, y_train)