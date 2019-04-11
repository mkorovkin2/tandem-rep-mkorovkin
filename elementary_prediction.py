import numpy as np
import pandas as pd
from scipy.stats import linregress

def show_hg_file2():
    file_path = "/Users/mkorovkin/Desktop/marzd/HG007counts.csv"

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
    df2.drop(['chrom', 'start', 'end'], axis=1)
    df2 = df2.loc[df2['indistinguishable'] == "S"]
    df2['TRID'] = df2['TRID'].apply(lambda x: str(x) + "_S")

    df = pd.merge(df, df2, on='TRID')

    x = df['array_length']

    return x, df['ratio'], df['pattern_size']

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return (s, i), (s1, i1)

def test_for_read_length(read_length):
    Alist2 = np.arange(read_length * 3, read_length * (40 - 7), read_length)

    mean_fit, std_fit = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length)
    test_x, test_y, test_patterns = show_hg_file2()

    in_range = 0
    less_count = 0
    greater_count = 0

    for i in range(len(test_x)):
        x = test_x[i]
        y = test_y[i]

        expected_y = mean_fit[0] * x + mean_fit[1]
        expected_std = std_fit[0] * x + mean_fit[1]
        expected_range = (expected_y - expected_std, expected_y + expected_std)

        if y > expected_range[1]:
            greater_count += 1
        elif y < expected_range[0]:
            less_count += 1
        else:
            in_range += 1

    total_range = greater_count + less_count + in_range

    in_range_perc = total_range / in_range
    less_range_perc = less_count / in_range
    greater_range_perc = greater_count / in_range

    print(in_range_perc, less_range_perc, greater_range_perc)

read_length_options = [100, 148, 250]

for read in read_length_options:
    test_for_read_length(read)