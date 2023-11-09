# Load C. diff. Dataset denoised and calculate mean counts over all Subjects, then save as csv file
import numpy as np
import pandas as pd
import pickle as pkl


data_path = "input_data/clv/bucci/"

Y = pkl.load(open(f"{data_path}Y_cdiff-denoised.pkl", "rb"))
U = pkl.load(open(f"{data_path}U_cdiff.pkl", "rb"))
T = pkl.load(open(f"{data_path}T_cdiff.pkl", "rb"))

# take mean over all subjects
Y = [np.sum(Y, axis=0)]
T = [T[0]]

# provide species names
names = [
    'Clostridium-hiranonis',
    'Clostridium-difficile',
    'Proteus-mirabilis',
    'Clostridium-scindens',
    'Ruminococcus-obeum',
    'Clostridium-ramosum',
    'Bacteroides-ovatus',
    'Akkermansia-muciniphila',
    'Parabacteroides-distasonis',
    'Bacteroides-fragilis',
    'Bacteroides-vulgatus',
    'Klebsiella-oxytoca',
    'Roseburia-hominis',
    'Escherichia-coli'
 ]

# save Bucci denoised mean data as csv file
df_data = pd.DataFrame(Y[0], columns=names)
df_data.insert(0, "Time", T[0])
df_data.to_csv(f'data/final_datasets/ts_bucci_subject_all_rel_counts_denoised.csv', index=False)