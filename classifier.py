import pandas as pd
import numpy as np
from scipy.signal import wiener
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import linregress

def load_hg_file(file_name):
    file_path = "/Users/mkorovkin/Desktop/marzd/" + file_name + ".csv"

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

    return df

def lin(dataset_x, dataset_y):
    # First definition of loss:
    b1h, b0h, _, _, _ = linregress(dataset_x, dataset_y)

    return b1h, b0h

def loss_of(dataset_eq, equation):
    b0 = equation[1]
    b1 = equation[0]

    b0h = dataset_eq[1]
    b1h = dataset_eq[0]

    return np.sinh(np.abs((b0 - b0h) / b0) + np.abs((b1 - b1h) / b1))

def loss_of_chi_sq(dataset_eq, equation):
    b0 = equation[1]
    b1 = equation[0]

    b0h = dataset_eq[1]
    b1h = dataset_eq[0]

    return np.abs((b0 - b0h) / b0) ** 2 + np.abs((b1 - b1h) / b1) ** 2

def hash_func(i, j):
    return 5 ** i + 5 ** j

def combinations_of_2(elements):
    set = {}

    for i in elements.keys():
        for j in elements.keys():
            if i != j:
                z = hash_func(i, j)
                if not (z in set):
                    set[z] = { i: elements[i], j: elements[j] }

    return set

def iterative_loss(dataset_x, dataset_y, equation_set, use_tuple=False):
    b1h, b0h = lin(dataset_x, dataset_y)

    deq = (b1h, b0h)

    eq_comb = combinations_of_2(equation_set)

    max_loss_origin = -1
    max_loss_origin_tuple = -1
    max_loss = -1

    for comb_key in eq_comb.keys():
        subset = eq_comb[comb_key]

        eqs = list()
        skey_list = list()
        for skey in subset.keys():
            eqs.append(subset[skey])
            skey_list.append(skey)

        loss1 = loss_of(deq, eqs[0])
        loss2 = loss_of(deq, eqs[1])
        nloss = loss1 + loss2
        if nloss > max_loss:
            max_loss = nloss
            if use_tuple:
                max_loss_origin_tuple = skey_list
            else:
                if loss1 > loss2:
                    max_loss_origin = skey_list[0]
                else:
                    max_loss_origin = skey_list[1]

    if use_tuple:
        del equation_set[max_loss_origin_tuple[0]]
        del equation_set[max_loss_origin_tuple[1]]
    else:
        del equation_set[max_loss_origin]

    return equation_set

def selection_loss(y_dataset_labels, x_selection, eq):
    predicted_values = x_selection * eq[0] + eq[1]

    loss = (np.abs(y_dataset_labels - predicted_values) ** 2).sum()

    return loss

def get_x_y_selection(dataset_x, dataset_y, selections):
    x_selection = np.array(np.random.uniform(0, len(dataset_x), selections), dtype=np.int32)

    y_dataset_labels = list()
    for x in x_selection:
        y_dataset_labels.append(dataset_y[x])
    y_dataset_labels = np.array(y_dataset_labels)

    return x_selection, y_dataset_labels

def iterative_loss_with_selection(dataset_x, dataset_y, equation_set, selections):
    x_selection, y_dataset_labels = get_x_y_selection(dataset_x, dataset_y, selections)

    max_loss = -1
    max_loss_key = -1

    for key in equation_set.keys():
        eq = equation_set[key]
        loss = selection_loss(y_dataset_labels, x_selection, eq)

        if loss > max_loss:
            max_loss = loss
            max_loss_key = key

    del equation_set[max_loss_key]

    return equation_set

def iterative_forwards_loss_with_selection(dataset_x, dataset_y, equation_set, selections):
    x_selection, y_dataset_labels = get_x_y_selection(dataset_x, dataset_y, selections)

    min_loss = 1000000
    min_loss_key = -1

    for key in equation_set.keys():
        eq = equation_set[key]
        loss = selection_loss(y_dataset_labels, x_selection, eq)

        if loss < min_loss:
            min_loss = loss
            min_loss_key = key

    return min_loss_key

def match_read_length(dataset_x, dataset_y):
    path = "/Users/mkorovkin/Desktop/marzd/output40_new.csv"
    mean_arrays = pd.read_csv(path)

    dataset_eq = linregress(dataset_x, dataset_y)

    keys_inspection = ["mA100", "mA148", "mA250"]

    result_dict = { 100: 0 }
    min_loss = 100

    for key in keys_inspection:
        dataf = np.array(mean_arrays[key])
        R = np.int32(key[2:])
        subAlist = np.arange(R * 3, R * 40, R)

        b1, b0, _, _, _ = linregress(subAlist, dataf)
        eq = (b1, b0)

        h_loss = np.int32(loss_of(dataset_eq, eq))

        result_dict[h_loss] = R

        if h_loss < min_loss:
            min_loss = h_loss

    return result_dict[min_loss]

def match_read_length_with_selection(dataset_x, dataset_y, selections):
    path = "/Users/mkorovkin/Desktop/marzd/output40_new.csv"
    mean_arrays = pd.read_csv(path)

    dataset_eq = linregress(dataset_x, dataset_y)

    keys_inspection = ["mA100", "mA148", "mA250"]

    result_dict = { 100: 0 }
    min_loss = 100

    x_selection, y_dataset_labels = get_x_y_selection(dataset_x, dataset_y, selections)

    for key in keys_inspection:
        dataf = np.array(mean_arrays[key])
        R = np.int32(key[2:])
        subAlist = np.arange(R * 3, R * 40, R)

        b1, b0, _, _, _ = linregress(subAlist, dataf)
        eq = (b1, b0)

        h_loss = np.int32(selection_loss(x_selection, y_dataset_labels, eq))

        result_dict[h_loss] = R

        if h_loss < min_loss:
            min_loss = h_loss

    return result_dict[min_loss]

def contained_in_area(y, s, yr):
    return yr > (y - s) and yr < (y + s)

def match_point_as_nearest_neighbor(point, total_set):
    path = "/Users/mkorovkin/Desktop/marzd/output40_new.csv"
    mean_arrays = pd.read_csv(path)

    significant_differences_labels = list()

    smallest_loss = 10000
    smallest_loss_label = "inconclusive"

    for (equation, std_eq, lab, _) in total_set:
        y_label = point['x'] * equation[0] + equation[1]
        y_std = point['x'] * std_eq[0] + std_eq[1]

        if not contained_in_area(y_label, y_std, point['y']):
            significant_differences_labels.append(lab)

        raw_loss = (y_label - point['y']) ** 2 / np.log(point['x'] / (point['i']) + 1)

        if smallest_loss > raw_loss and (not (lab in significant_differences_labels)):
            smallest_loss = raw_loss
            smallest_loss_label = lab

    return significant_differences_labels, smallest_loss_label

def scaled_homozygous_loss(y_label, point, operation):
    return (y_label - point['y'] * (point['y'] + operation * point['p']) / point['y']) ** 2

def match_point_assume_pattern_nearest_neighbor(point, total_set):
    path = "/Users/mkorovkin/Desktop/marzd/output40_new.csv"
    mean_arrays = pd.read_csv(path)

    significant_differences_labels = list()

    smallest_loss = 10000
    smallest_loss_label = "inconclusive"

    for (equation, std_eq, lab, operation) in total_set:
        y_label = point['x'] * equation[0] + equation[1]
        y_std = point['x'] * std_eq[0] + std_eq[1]

        if not contained_in_area(y_label, y_std, point['y']):
            significant_differences_labels.append(lab)

        raw_loss = (y_label - point['y']) ** 2 * scaled_homozygous_loss(y_label, point, operation)

        if smallest_loss > raw_loss and (not (lab in significant_differences_labels)):
            smallest_loss = raw_loss
            smallest_loss_label = lab

    return significant_differences_labels, smallest_loss_label

def run_script(file_index):
    # Load file
    file_list = ["HG001counts",
                 "HG002counts",
                 "HG006counts",
                 "HG007counts"
                 ]
    df = load_hg_file(file_list[file_index])

    estimated_read = match_read_length_with_selection(df['array_length'], df['ratio'], 2000)
    print("R=" + str(estimated_read))

    print("---")
    print("Loss computation results (not very relevant):")

    if estimated_read == 100:
        eq00 = (0.008947817221045604, -0.9989622172344781)
        eq01 = (0.010125161592235813, -1.1304046142390156)
        eq0n1 = (0.007770472849855392, -0.8675198202299388)
        eq11 = (0.011302505963426022, -1.2618470112435425)
        eqn11 = (0.006593128478665182, -0.7360774232254066)
        st = (0.0013359142445383099, -0.026017551104345227)
    elif estimated_read == 148:
        eq00 = (0.006008077491641273, -0.98365355261177)
        eq01 = (0.006798614003699335, -1.1130816516396251)
        eq0n1 = (0.005217540979583209, -0.8542254535838989)
        eq11 = (0.007589150515757398, -1.2425097506674945)
        eqn11 = (0.004427004467525148, -0.7247973545560331)
        st = (0.0008079342644219522, -0.11799680124040313)
    else:
        eq00 = (0.003561945691532275, -0.9703238049866592)
        eq01 = (0.004030622756207573, -1.0979979898533152)
        eq0n1 = (0.0030932686268569745, -0.8426496201199853)
        eq11 = (0.004499299820882874, -1.2256721747199926)
        eqn11 = (0.0026245915621816752, -0.7149754352533204)
        st = (0.0003861314031995404, -0.10155638596216687)

    b1h, b0h = lin(df['array_length'], df['ratio'])
    dataset_eq = (b1h, b0h)

    for (eq, lab) in [(eq00, "0/0"), (eq01, "0/1"), (eq11, "1/1"), (eq0n1, "0/-1"), (eqn11, "-1/-1")]:
        ba = [250, 1000]
        # plt.plot(ba, [eq[1] + x * eq[0] for x in ba])

        print(lab, "=", np.round(loss_of(dataset_eq, eq), decimals=4))
    # plt.show()

    print("---")

    equation_set = {0: eq00, 1: eq01, 2: eq11, 3: eq0n1, 4: eqn11}
    equation_set_labels = {0: "0/0", 1: "0/1", 2: "1/1", 3: "0/-1", 4: "-1/-1"}
    while len(equation_set.keys()) > 1:
        equation_set = iterative_loss_with_selection(df['array_length'], df['ratio'], equation_set, 1500)
    for k in equation_set.keys():
        print("Estimated genotype:", equation_set_labels[k])

    correct_genotype_estimated = equation_set_labels[iterative_forwards_loss_with_selection(df['array_length'], df['ratio'], equation_set, 1500)]
    print("Estimated genotype (based on forward selection:", correct_genotype_estimated)

    print("---")

    total_set = [(eq00, st, "0/0", 0),
                 (eq01, st, "0/1", 0.5),
                 (eq11, st, "1/1", 1),
                 (eq0n1, st, "0/-1", -0.5),
                 (eqn11, st, "-1/-1", -1)]

    # Select a few random points from the HG001 genome
    random_rows = np.array(np.random.uniform(0, len(df['TRID']), 1500), dtype=np.int32)

    result_keys_df = {}
    for row in random_rows:
        point = {
            "x": df.iloc[row]['array_length'],
            "y": df.iloc[row]['ratio'],
            "c": df.iloc[row]['copy_number'],
            "p": df.iloc[row]['pattern_size'],
            "i": df.iloc[row]['inside_array'],
            "f": df.iloc[row]['right_flank'] + df.iloc[row]['left_flank'],
        }
        _, genotype = match_point_as_nearest_neighbor(point, total_set)

        if not (genotype in result_keys_df):
            result_keys_df[genotype] = {}
            result_keys_df[genotype]['value'] = 1
            result_keys_df[genotype]['patterns'] = {}
            result_keys_df[genotype]['patterns'][np.int32(point['p'])] = 1
            result_keys_df[genotype]['xs'] = {}
            result_keys_df[genotype]['xs'][np.int32(point['x'])] = 1
        else:
            result_keys_df[genotype]['value'] = result_keys_df[genotype]['value'] + 1
            ps_temp = np.int32(point['p'])
            if ps_temp in result_keys_df[genotype]['patterns']:
                result_keys_df[genotype]['patterns'][ps_temp] = result_keys_df[genotype]['patterns'][ps_temp] + 1
            else:
                result_keys_df[genotype]['patterns'][ps_temp] = 1

            x_temp = np.int32(point['x'])
            if x_temp in result_keys_df[genotype]['xs']:
                result_keys_df[genotype]['xs'][x_temp] = result_keys_df[genotype]['xs'][x_temp] + 1
            else:
                result_keys_df[genotype]['xs'][x_temp] = 1

    del result_keys_df["inconclusive"]
    perc_total_sum_count = 0

    for genotype in result_keys_df.keys():
        genotype_df = result_keys_df[genotype]

        print(genotype, "-> total values:", result_keys_df[genotype]['value'])
        perc_total_sum_count += result_keys_df[genotype]['value']

        pat_list = list()
        for k in result_keys_df[genotype]['patterns'].keys():
            pat_list.append(k)

        x_list = list()
        for k in result_keys_df[genotype]['xs'].keys():
            x_list.append(k)

        p_mean = np.int32(np.array(pat_list).mean())
        print("-----> P-mean:", p_mean)

        a_mean = np.int32(np.array(x_list).mean())
        print("-----> A-mean:", a_mean)

        print("-----> Result:", np.round(a_mean / p_mean, decimals=2))

    perc_correct = result_keys_df[correct_genotype_estimated]['value'] / perc_total_sum_count
    print("Correct percentage of", correct_genotype_estimated, "=", str(np.round(perc_correct * 100.0, decimals=2)) + "%")

    ########################################################################################################################

    exit(0)

    def get_two_random_rows(df):
        rand_nums = np.random.randint(0, len(df['TRID']), 2)

        return df.iloc[rand_nums[0]], df.iloc[rand_nums[1]]

    # Select two random rows from the dataframe
    row1, row2 = get_two_random_rows(df)

    def are_they_different(row1, row2):
        pass # row1['ratio'], row2['ratio']

    # Can we tell that they're different just based on ratio?
    different_result = are_they_different(row1, row2)
    print(different_result)

    # ['TRID', 'chromosome', 'start_x', 'end_x', 'array_length',
    # 'inside_array', 'left_flank', 'right_flank', 'Unnamed: 8',
    # 'ratio', 'indistinguishable', 'chrom', 'start_y', 'end_y',
    # 'pattern_size', 'array_length_refset', 'copy_number']

    print(list(df))
    print(df.iloc[11])