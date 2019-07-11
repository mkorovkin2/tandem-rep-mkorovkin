import pandas as pd
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])


    print(len(Alist))
    print(len(means))
    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return np.array([(s * A + i) for A in Alist]), np.array([(s1 * A + i1) for A in Alist])

file_path = "/Users/mkorovkin/Desktop/marzd/HG007counts.csv"

df = pd.read_csv(file_path)

df = df.loc[df['left_flank'] > 0]
df = df.loc[df['right_flank'] > 0]
df = df.loc[df['array_length'] > 0]
df = df.loc[df['inside_array'] > 0]
df = df.loc[df['array_length'] < 5000]

df['ratio'] = df['inside_array'] / (df['left_flank'] + df['right_flank'])

x = df['array_length']

plt.scatter(x, df['ratio'])
read_length = 100
Alist = np.arange(read_length * 3, read_length * (40 - 7), read_length)
mean, std = get_mean_std("/Users/mkorovkin/Desktop/marzd/output30.csv", Alist, read_length)
plt.plot(Alist, mean)
plt.plot(Alist, mean + std * 1.94)
plt.plot(Alist, mean - std * 1.94)
plt.show()