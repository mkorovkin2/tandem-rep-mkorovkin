import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

def pandas_analyze(path, expected_value=0, expected_values=[], genotype_filter="", printE=0, title=""):
    dataframe = pd.read_csv(path)
    #['TRID', 'chromosome', 'start', 'end', 'array_length', 'inside_array', 'left_flank', 'right_flank', 'ratio_array']
    # seems like inside_array/array_length = ratio_array
    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]
    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]

    dataframe = dataframe.loc[dataframe['array_length'] < 10000]

    inside_array = np.array(dataframe.iloc[:]["inside_array"])
    right_flank = np.array(dataframe.iloc[:]["right_flank"])
    left_flank = np.array(dataframe.iloc[:]["left_flank"])
    array_length = np.array(dataframe.iloc[:]["array_length"])

    ratio_list = inside_array / (right_flank + left_flank)

    plt.scatter(array_length, ratio_list)
    plt.xlabel("Total length of analyzed sequences")
    plt.ylabel("Ratio X/Y")

    slope, intercept, _, _, _ =linregress(array_length, ratio_list)
    print("- - - - -")
    print(slope * 100, "per 100 |", intercept, "|", intercept + slope * 100)
    print(slope * 250, "per 250 |", intercept, "|", intercept + slope * 250)
    if expected_value > 0:
        print("Expected value:", slope * expected_value + intercept)
    elif len(expected_values) > 0:
        for e in expected_values:
            print("Expected value:", slope * e + intercept, "for", e)
    if len(title) > 0:
        plt.title(
            "READ_LENGTH = 250\nFor  an array length of 1750, [expected ratio = " + str(np.round(slope * printE + intercept, decimals=2)) + "] | " + title)
    else:
        plt.title(
            "For  an array length of 700, [expected ratio = " + str(np.round(slope * printE + intercept, decimals=2)) + "]")
    plt.show()
    print()

def analyze(path):
    content = list()
    with open(path) as f:
        for line in f:
            content.append(line.strip().split())

    flank1_read = {}
    flank1_bp = {}
    flank2_read = {}
    flank2_bp = {}
    array_read = {}
    array_bp = {}
    unique_keys = list()
    remove_keys = list()

    for item in content:
        try:
            print(item)
            break
        except:
            print("Error")

    print(len(unique_keys))
    sum_ratio = 0
    count = 0
    ratio_dict = {}
    dlist = list()
    upper = 1500
    lower = 0
    for key in unique_keys:
        if (flank1_read[key] + flank2_read[key]) > 0 and array_read[key] < upper and array_read[key] > lower:
            ratio = array_read[key] / (flank1_read[key] + flank2_read[key])
            if ratio > 10:
                pass
            dlist.append(ratio)
            ratio_dict[ratio] = (array_bp[key] + flank1_bp[key] + flank2_bp[key])
            sum_ratio += ratio
            count += 1
    return ratio_dict

#ratio_dict = pandas_analyze("/Users/mkorovkin/Desktop/bp_analysis-228/simulations - 100bp.csv")
#ratio_dict = pandas_analyze("/Users/mkorovkin/Desktop/bp_analysis-228/simulations - 148bp.csv")
#ratio_dict = pandas_analyze("/Users/mkorovkin/Desktop/bp_analysis-228/simulations - 250bp.csv", expected_values=[500, 1000, 1500, 2000])
ratio_dict1 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/sim1_250bp.csv", expected_values=[500, 700, 1000], genotype_filter='0/0', printE=1750, title="0/0")
ratio_dict2 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/sim1_250bp.csv", expected_values=[500, 700, 1000], genotype_filter='0/1', printE=1750, title="0/1")
ratio_dict3 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/sim1_250bp.csv", expected_values=[500, 700, 1000], genotype_filter='0/-1', printE=1750, title="0/-1")
ratio_dict4 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/sim1_250bp.csv", expected_values=[500, 700, 1000], genotype_filter='43466', printE=1750, title="1/1")