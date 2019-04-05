import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt

def find_variance(kwargs, A_v, R_v, F_v):
    to_find = "variance"
    followup_variance = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip()
    if "in" in followup_variance:
        to_find = "in"
        followup_in = followup_variance[(followup_variance.index(to_find) + len(to_find) + 1):].strip()
        if "A" in followup_in:
            to_find = "A"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            A_v = (True, np.float32(followup[0]))
        elif "R" in followup_in:
            to_find = "R"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            R_v = (True, np.float32(followup[0]))
        elif "F" in followup_in:
            to_find = "F"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            F_v = (True, np.float32(followup[0]))
    if "variance" in followup_variance:
        return find_variance(followup_variance, A_v, R_v, F_v)
    else:
        return A_v, R_v, F_v

def sim(kwargs):
    loss = (False, 0, 0)
    gain = (False, 0, 0)
    normal = (False, 0, 0)
    A_v = (False, 0)
    R_v = (False, 0)
    F_v = (False, 0)

    if len(kwargs) > 0:
        if "loss" in kwargs:
            to_find = "loss"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            loss = (True, np.float32(followup[0]), np.float32(followup[1]))
        if "gain" in kwargs:
            to_find = "gain"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            gain = (True, np.float32(followup[0]), np.float32(followup[1]))
        if "normal" in kwargs:
            to_find = "normal"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            normal = (True, np.float32(followup[0]), np.float32(followup[1]))
        if "variance" in kwargs:
            A_v, R_v, F_v = find_variance(kwargs, A_v, R_v, F_v)

    print(loss)
    print(gain)
    print(normal)
    print(A_v)
    print(R_v)
    print(F_v)

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)

    return [(s * A + i) for A in Alist]

read_length = 100
Alist = np.arange(read_length * 3, read_length * 40, read_length)
gs = get_mean_std("/Users/mkorovkin/Desktop/marzd/output30.csv", Alist, read_length)

print(gs)

plt.scatter(Alist, gs)
plt.show()