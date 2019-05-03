import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt

def show_hg_file2(file_path):
    df = pd.read_csv(file_path)

    df = df.loc[df['left_flank'] > 0]
    df = df.loc[df['right_flank'] > 0]
    df = df.loc[df['array_length'] > 0]
    df = df.loc[df['inside_array'] > 0]
    df = df.loc[df['array_length'] < 5000]
    df = df.loc[df['TRID'].str.contains("S")]

    df['ratio'] = df['inside_array'] / (df['left_flank'] + df['right_flank'])
    df = df.loc[df['ratio'] < 200]

    file_path2 = "/Users/mkorovkin/Desktop/marzd/refset_full.csv"

    df2 = pd.read_csv(file_path2)
    # print(list(df2))
    df2.drop(['chrom', 'start', 'end'], axis=1)
    df2 = df2.loc[df2['indistinguishable'] == "S"]
    df2['TRID'] = df2['TRID'].apply(lambda x: str(x) + "_S")

    df = pd.merge(df, df2, on='TRID')

    x = df['array_length']

    return x, df['ratio'], df['pattern_size'], df['copy_number']

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return (s, i), (s1, i1)

def test_for_read_length(read_length, file_path):
    Alist2 = np.arange(read_length * 3, read_length * (40 - 7), read_length)

    # -> output40_coverage355
    mean_fit, std_fit = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length)
    test_x, test_y, test_patterns, test_copy_number = show_hg_file2(file_path)

    test_patterns = np.log(test_patterns)

    in_range = 0
    less_count = 0
    greater_count = 0

    for i in range(len(test_x)):
        x = test_x[i]
        y = test_y[i]

        expected_y = mean_fit[0] * x + mean_fit[1]
        expected_std = std_fit[0] * x + mean_fit[1]
        expected_range = (expected_y - expected_std, expected_y + expected_std)

        scale = test_patterns[i]

        if y > expected_range[1]:
            greater_count += 1 / scale
        elif y < expected_range[0]:
            less_count += 1 / scale
        else:
            in_range += 1 / scale

    total_range = greater_count + less_count + in_range

    in_range_perc = in_range / total_range
    less_range_perc = less_count / total_range
    greater_range_perc = greater_count / total_range

    ratio_basic = max(less_range_perc, greater_range_perc) / min(less_range_perc, greater_range_perc)
    print(str(np.round(less_range_perc, decimals=4)) + "<--" + str(np.round(in_range_perc, decimals=4)) + "-->" + str(np.round(greater_range_perc, decimals=4)), ratio_basic)

    # plt.scatter(test_patterns, test_copy_number)
    # plt.show()

def look_for():
    read_length_options = [100, 148, 250]

    file_paths = [
        "/Users/mkorovkin/Desktop/marzd/HG001counts.csv",
        "/Users/mkorovkin/Desktop/marzd/HG002counts.csv",
        "/Users/mkorovkin/Desktop/marzd/HG007counts.csv"

    ]
    labels = [
        148,
        250,
        148
    ]

    for file_index in range(len(file_paths)):
        for read in read_length_options:
            test_for_read_length(read, file_paths[file_index])
        print("Label: " + str(labels[file_index]))
        print("---")
        break

#look_for()
#test_for_read_length(148, "/Users/mkorovkin/Desktop/marzd/HG007counts.csv")
# test_x, test_y, test_patterns, test_copy_number = show_hg_file2("/Users/mkorovkin/Desktop/marzd/HG007counts.csv")



#read length/coverage vs. standard deviation
# Title graphs with coverage and read length/ indicate that it/s 0/0 or something
# Graph the flank coverage (graph the errors against the flank)
# Find expected flanks, and compute a regression line for expected flanks

# Aggregate graphs and different versinos of what we're lookign at
    # Points are above the line
    # Flanks seem to have an effect on -> distribution of the flank values for each array length
    # Distribution of inside array length for each array length
    # -> See which one is causing the problem
        # Use average value for flanks (across the whole genome)

# Put int MD format/notebook
# Smooth standard deviation
# INtervals for standard deviation for a number of points versus X axis
# Make the flanks a static number; that's why small patterns seem to jump above
    # Getting a flanks distributions (use mean value)

# Exclude certain lower pattern lengths

# Put into a python notebook