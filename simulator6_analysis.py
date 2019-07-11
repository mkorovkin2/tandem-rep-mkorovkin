import pandas as pd
import numpy as np
from scipy.stats import linregress

df = pd.read_csv("/Users/mkorovkin/Desktop/output_statistics_00_new.csv")

means = [np.array(df[['mA100']]).transpose(),
         np.array(df[['mA148']]).transpose(),
         np.array(df[['mA250']]).transpose()]
stds = [np.array(df[['sA100']]).transpose(),
        np.array(df[['sA148']]).transpose(),
        np.array(df[['sA250']]).transpose()]

means = np.array(means)
stds = np.array(stds)

means_plus_stds = 2 * stds

for i in range(len(means)):
    # This value is how different ratios must be to determine a statistically significant difference
    print(means_plus_stds)