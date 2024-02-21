import os
import numpy as np
import pickle
import pandas as pd


sample_size = 90

disc_muts_full = np.linspace(0,1,sample_size+1)[1:]

with open(f'./disc_muts_full{sample_size}.pkl', 'wb') as handle:
    pickle.dump(disc_muts_full, handle)

with open(f'./times.pkl', 'rb') as handle:
    times = pickle.load(handle)
with open(f'./thetas.pkl', 'rb') as handle:
    thetas = pickle.load(handle)
with open(f'./disc_muts_full{sample_size}.pkl', 'rb') as handle:
    disc_muts_full = pickle.load(handle)
    
    
loger2 = np.zeros((10, times.shape[0], disc_muts_full.shape[0], thetas.shape[0]))
for mut in [0,1,2,3,4,5,6,7,8,9]:
    for time in range(0,800):
        
        loger2[mut, time] = np.load(f'./Likelihoods/{sample_size}M={mut}/{time}.pkl', allow_pickle=True)

loger2 = np.swapaxes(loger2, 1, 2)

with open(f'./logLikelySingleCell{sample_size}.pkl', 'wb') as handle:
    pickle.dump(loger2, handle)