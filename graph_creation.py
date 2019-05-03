import pandas as pd
import numpy as np
from scipy.stats import linregress
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as filter

def scale_to_range(list_, min_, max_):
    list_ = np.log(list_)
    list_ = np.array(list_, dtype=np.float32)
    min_sp = list_.min()
    max_sp = list_.max()
    list_ = list_ - min_sp
    list_ = list_ / max_sp
    list2 = (1 - list_) * min_ / max_
    list_ = (list_ + list2) * max_
    return np.array(np.round(list_, decimals=0), dtype=np.int32)

def upper_percentile_ps(df, perc):
    perc = np.int32(np.max(np.floor(perc * 100), 0))
    breakp = np.percentile(np.array(df['pattern_size'], dtype=np.float32), perc)
    return df.loc[df['pattern_size'] > breakp]

def upper_percentile_al(df, perc):
    perc = np.int32(np.max(np.floor(perc * 100), 0))
    breakp = np.percentile(np.array(df['array_length'], dtype=np.float32), perc)
    return df.loc[df['array_length'] > breakp]
# print(scale_to_range([1, 2, 3, 4, 5, 6, 7, 8], 0, 39-15))

def invert_colors(array):
    minx = array.min()
    maxx = array.max()

    return array * -1 + maxx + minx

def show_array_length_frequency_plot(data0, infix):
    plt.hist(data0, bins=40)
    plt.xlabel("Array lengths in dataset")
    plt.ylabel("Frequency of array lengths")
    plt.title("Array Length Frequency of " + infix + " Data")
    plt.show()

def pandas_analyze(path, READ_LENGTH, genotype_filter="", printE=0, title="", force_exclude=0, force_stop=0,
                   drop_lowest_percent_pattern_length=0.0, drop_lowest_percent_array_length=0.0):
    dataframe = pd.read_csv(path)
    if force_stop != 0:
        dataframe = dataframe.iloc[:force_stop]

    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    else:
        dataframe = dataframe.loc[dataframe['array_length'] < 15000]

    if drop_lowest_percent_array_length > 0 and drop_lowest_percent_array_length < 1:
        dataframe = upper_percentile_al(dataframe, drop_lowest_percent_array_length)

    if drop_lowest_percent_pattern_length > 0 and drop_lowest_percent_pattern_length < 1:
        dataframe = upper_percentile_ps(dataframe, drop_lowest_percent_pattern_length)

    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]
    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]
    dataframe['ratio'] = dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank'])

    scaled_patterns = dataframe['pattern_size']
    scaled_patterns = scale_to_range(scaled_patterns, 0, 39 - 15)

    return dataframe, scaled_patterns, np.array(dataframe['array_length'])

def score_df(csv_name, minx, maxx, xs, force_exclude=0, static=(False, 0)):
    dataframe = pd.read_csv(csv_name)
    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    dataframe = dataframe.loc[dataframe['right_flank'] != 0]
    dataframe = dataframe.loc[dataframe['genotype'] == '0/0']
    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    array_len_list = np.array(dataframe['array_length'])
    ratio_list = np.array(dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank']))

    slope, rint, _, _, _ = linregress(array_len_list, ratio_list)

    if not static[0]:
        print("Regression line equation: y = (" + str(np.round(rint, decimals=4)) + ") + (" + str(
            np.round(slope, decimals=4)) + ")x")
        return [minx, maxx], [slope * minx + rint, slope * maxx + rint], np.array([(slope * x + rint) for x in xs])
    else:
        if static[1] == 100:
            return [minx, maxx], [0.0118914 * minx + -2.6604506, 0.0118914 * maxx + -2.6604506], np.array([(0.0118914 * x + -2.6604506) for x in xs])
            #return [minx, maxx], [0.0049349 * minx + -1.988381, 0.0049349 * maxx + -1.988381], np.array([(0.0049349 * x + -1.988381) for x in xs])

        elif static[1] == 148:
            return [minx, maxx], [0.008597 * minx + -2.7313725, 0.008597 * maxx + -2.7313725], np.array([(0.008597 * x + -2.7313725) for x in xs])
            #return [minx, maxx], [0.0084543 * minx + -0.8394286, 0.0084543 * maxx + -0.8394286], np.array([(0.0084543 * x + -0.8394286) for x in xs])

        elif static[1] == 250:
            return [minx, maxx], [0.0053299 * minx + -2.4744479, 0.0053299 * maxx + -2.4744479], np.array([(0.0053299 * x + -2.4744479) for x in xs])
            #return [minx, maxx], [0.0065058 * minx + -1.2928571, 0.0065058 * maxx + -1.2928571], np.array([(0.0065058 * x + -1.2928571) for x in xs])

        else:
            return
    #return [minx, maxx], [slope * minx + rint, slope * maxx + rint], np.array([(slope * x + rint) for x in xs])

def score_df_nofilter(csv_name, minx, maxx, xs, force_exclude=0):
    dataframe = pd.read_csv(csv_name)
    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    dataframe = dataframe.loc[dataframe['right_flank'] != 0]
    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    array_len_list = np.array(dataframe['array_length'])
    ratio_list = np.array(dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank']))

    slope, rint, _, _, _ = linregress(array_len_list, ratio_list)

    print("Regression line equation: y = (" + str(np.round(rint, decimals=4)) + ") + (" + str(np.round(slope, decimals=4)) + ")x")
    return [minx, maxx], [slope * minx + rint, slope * maxx + rint], np.array([(slope * x + rint) for x in xs])

def score_df_for_std_from_file(bp, coverage, regress_std=False):
    Alist = np.arange(bp * 3, bp * 50, bp)
    temp_df = pd.read_csv("/Users/mkorovkin/Desktop/marzd/std_" + str(bp) + ".csv")
    std_list = np.array(temp_df[[str(coverage) if (coverage == 70 or coverage == 355) else "def"]],
                        dtype=np.float32)
    if not regress_std:
        if bp == 100 or bp == 148 or bp == 250:
            return Alist, std_list.transpose()[0] # [minx, maxx], np.array([minx * eq[0] + eq[1], maxx * eq[0] + eq[1]])
    else:
        if bp == 100 or bp == 148 or bp == 250:
            split_point = 3
            std_list2 = list(std_list[:split_point])
            std_reg = std_list[split_point:]
            std_reg = std_reg.reshape(len(std_reg),)
            mslope, mintercept, _, _, _ = linregress(Alist[split_point:], std_reg)
            for A in Alist[split_point:]:
                std_list2.append(mslope * A + mintercept)
            return Alist, np.array(std_list2, dtype=np.float32) # [minx, maxx], np.array([minx * eq[0] + eq[1], maxx * eq[0] + eq[1]])
    return
    #(0.00531795744680851, -1.5549748226950353) <- for C=70
    #(0.0022562765957446804, -0.8927053191489343) <- for C=355
    #(0.0037871170212765964, -1.2238400709219857) <- for MEAN

def pandas_analyze_new(path, READ_LENGTH, genotype_filter="", printE=0, title="",
                       force_exclude=0, force_stop=0, drop_lowest_percent_pattern_length=0.0,
                       drop_lowest_percent_array_length=0.0, copy_number_exclusion=0):
    dataframe = pd.read_csv(path)
    if force_stop != 0:
        dataframe = dataframe.iloc[:force_stop]

    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    else:
        dataframe = dataframe.loc[dataframe['array_length'] < 15000]

    if drop_lowest_percent_array_length > 0 and drop_lowest_percent_array_length < 1:
        dataframe = upper_percentile_al(dataframe, drop_lowest_percent_array_length)

    if drop_lowest_percent_pattern_length > 0 and drop_lowest_percent_pattern_length < 1:
        dataframe = upper_percentile_ps(dataframe, drop_lowest_percent_pattern_length)

    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]

    if copy_number_exclusion > 0:
        dataframe = dataframe.loc[dataframe['copy_number'] < copy_number_exclusion]

    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]
    dataframe['ratio'] = dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank'])

    scaled_patterns = dataframe['pattern_size']
    scaled_patterns = scale_to_range(scaled_patterns, 0, 39 - 15)

    s, i, _, _, stderr = linregress(dataframe['array_length'], dataframe['ratio'])
    Alist = np.arange(read_length * 3, read_length * 50, read_length)
    regression_array = np.array([(s * x + i) for x in Alist])

    error_bars = error_bars_around_dataset(dataframe, read_length * 3, read_length * 50, read_length)

    return dataframe, np.log(5 + scaled_patterns), error_bars, regression_array

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return np.array([(s * A + i) for A in Alist]), np.array([(s1 * A + i1) for A in Alist])

def error_bars_around_dataset(df, low, high, read_length):
    Alist3 = np.arange(low, high + read_length, read_length)

    std_list = list()
    for i in range(len(Alist3) - 1):
        minv = Alist3[i]
        maxv = Alist3[i + 1]

        df2 = df.loc[df['array_length'] <= maxv]
        df2 = df2.loc[df2['array_length'] >= minv]

        if len(df2) > 4:
            temp_std = np.array(df2['ratio']).std() * 1.96
            std_list.append(temp_std)
        else:
            std_list.append(0)

    return np.array(std_list)

local_path = "/Users/mkorovkin/Desktop/graphs_marzd/"

read_length = 250
Alist2 = np.arange(read_length * 3, read_length * (40 - 7), read_length)
Alist = np.arange(read_length * 3, read_length * 50, read_length)
coverage = 100#np.floor((355 - 70) / 2)

minax = 0
maxax = read_length * 50

csv_filename = ""
if read_length == 100:
    csv_filename = "sim3_100bp.csv"
elif read_length == 148:
    csv_filename = "sim2_148bp.csv"
elif read_length == 250:
    csv_filename = "sim1_250bp.csv"
else:
    exit(0)

excludeperc_ps = -0.7
excludeperc_arraylength = -0.9

alist_standard, stdlist_standard = score_df_for_std_from_file(read_length, 100)
_, _, ypoints = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, alist_standard, force_exclude=maxax, static=(True, read_length))

copy_number_exclude = 8 #12

df_00, color00, arraylength_data00, regression_arr00 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                          genotype_filter='0/0', force_exclude=maxax,
                                                                          drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                          drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                          copy_number_exclusion=copy_number_exclude)
df_0n1, color0n1, arraylength_data0n1, regression_arr0n1 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                              genotype_filter='0/-1', force_exclude=maxax,
                                                                              drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                              drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                                copy_number_exclusion=copy_number_exclude)
df_01, color01, arraylength_data01, regression_arr01 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                          genotype_filter='0/1', force_exclude=maxax,
                                                                          drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                          drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                          copy_number_exclusion=copy_number_exclude)
df_11, color11, arraylength_data11, regression_arr11 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                          genotype_filter='43466', force_exclude=maxax,
                                                                          drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                          drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                          copy_number_exclusion=copy_number_exclude)
df_n11, colorn11, arraylength_datan11, regression_arrn11 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                              genotype_filter='1', force_exclude=maxax,
                                                                              drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                              drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                                copy_number_exclusion=copy_number_exclude)
fitted_means_expected, fitted_stds = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length)

def show_std_plot():
    g = list()
    for read in [100, 148, 250]:
        fitted_means_expected, fitted_stds = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read)
        g.append(fitted_stds)

    trace0 = go.Scatter(
        x=Alist2,
        y=g[0],
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="Standard deviation for read=100"
    )
    trace1 = go.Scatter(
        x=Alist2,
        y=g[1],
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="Standard deviation for read=148"
    )
    trace2 = go.Scatter(
        x=Alist2,
        y=g[2],
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="Standard deviation for read=250"
    )
    layout = go.Layout(title="READ_LENGTHS of [100, 148, 250] | coverage=" + str(np.int32(coverage)) + " | array lengths from [" + str(read_length * 3) + ", " + str(read_length * 50) + "]",
                       xaxis=dict(ticks='', showticklabels=True,
                                  zeroline=True, title='Array length (base pairs)',
                                    titlefont=dict(
                                        family='Courier New, monospace',
                                        size=18,
                                        color='#7f7f7f'
                                    )),
                       yaxis=dict(ticks='', showticklabels=True,
                                  zeroline=True, title='Ratio (X/Y)',
                                    titlefont=dict(
                                        family='Courier New, monospace',
                                        size=18,
                                        color='#7f7f7f')),
                       showlegend=True, hovermode='closest', legend=dict(orientation='h'))
    py.plot(go.Figure([trace0, trace1, trace2], layout), filename='scatter-plot-with-colorscale')

def get_std(df, partitions):
    sorted_df = df.sort_values(['array_length'])

    interval = np.floor(len(sorted_df) / partitions)

    stda = list()
    stdx = list()

    for p in range(partitions):
        li = sorted_df['ratio'][np.int32(p * interval):np.int32((p + 1) * interval)]
        li_array = sorted_df['array_length'][np.int32(p * interval):np.int32((p + 1) * interval)]
        stda.append(np.std(li))
        stdx.append(max(li_array))

    li2 = sorted_df['ratio'][np.int32(partitions * interval):]
    li_array2 = sorted_df['array_length'][np.int32(partitions * interval):]
    stda.append(np.std(li2))
    stdx.append(max(li_array2))

    return np.array(stda), np.array(stdx)

def plot_stds_against_stds(df, bct, sed, appended):
    stda, stdx = get_std(df, sed)
    print(stdx)

    dfc = pd.DataFrame(data={"stda": stda, "stdx": stdx})
    # dfc.to_csv(local_path + "std_df_" + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".csv")

    plt.plot(Alist2, fitted_stds)
    plt.plot(stdx, filter_it(stda))
    plt.xlabel("Array length")
    plt.ylabel("Ratio standard deviation")
    plt.title("Array length vs. ratio standard deviation for READ_LENGTH=" + str(read_length) + "\nIntervals sampled: " + str(sed) +
              " | points per interval: " + str(np.int32(len(df['ratio']) / sed)) + "\n" + bct)
    # plt.show() #"cne" + str(copy_number_exclude) +
    # plt.savefig(local_path + "std_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def get_means(df, partitions):
    sorted_df = df.sort_values(['array_length'])

    interval = np.floor(len(sorted_df) / partitions)

    meansa = list()
    meansx = list()

    for p in range(partitions):
        li = sorted_df['ratio'][np.int32(p * interval):np.int32((p + 1) * interval)]
        li_array = sorted_df['array_length'][np.int32(p * interval):np.int32((p + 1) * interval)]
        meansa.append(np.mean(li))
        meansx.append(max(li_array))

    li2 = sorted_df['ratio'][np.int32(partitions * interval):]
    li_array2 = sorted_df['array_length'][np.int32(partitions * interval):]
    meansa.append(np.mean(li2))
    meansx.append(max(li_array2))

    return np.array(meansa), np.array(meansx)

def filter_it(npa):
    return npa # return filter(npa, window_length=3, polyorder=2)

def plot_means_against_means(df, bct, sed, appended):
    meansa, meansx = get_means(df, sed)

    plt.plot(Alist2, fitted_means_expected)
    plt.plot(meansx, filter_it(meansa))
    plt.xlabel("Array length")
    plt.ylabel("Ratio means")
    plt.title("Array length vs. ratio means for READ_LENGTH=" + str(read_length) + "\nIntervals sampled: " + str(sed) +
              " | points per interval: " + str(np.int32(len(df['ratio']) / sed)) + "\n" + bct)
    # plt.show()
    # plt.savefig(local_path + "mean_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def mean_and_std_figures():
    for (df, bct, appended) in [(df_n11, "Homozygous -1/-1", "n11"), (df_0n1, "Heterozygous 0/-1", "0n1"),
                                (df_00, "Homozygous 0/0", "00"), (df_01, "Heterozygous 0/1", "01"),
                                (df_11, "Homozygous 1/1", "11")]:
        sed = 50
        plot_stds_against_stds(df, bct, sed, appended)
        plot_means_against_means(df, bct, sed, appended)

def flank_figures():
    for (df, bct, appended) in [(df_n11, "Homozygous -1/-1", "n11"), (df_0n1, "Heterozygous 0/-1", "0n1"),
                                (df_00, "Homozygous 0/0", "00"), (df_01, "Heterozygous 0/1", "01"),
                                (df_11, "Homozygous 1/1", "11")]:
        flank_plotting(df, bct, appended, np.mean(df['right_flank'] + df['left_flank']), np.median(df['right_flank'] + df['left_flank']))

def flank_plotting(df, bct, appended, df_mean, df_median):
    plt.hist((df['right_flank'] + df['left_flank']), bins=20)
    plt.xlabel("Sum of flank counts")
    plt.ylabel("Frequency")
    plt.title("Histogram plot of flanks for READ_LENGTH=" + str(read_length) + "\n" + bct + "\nMean=" + str(np.round(df_mean, decimals=2)) + " | median=" + str(np.round(df_median, decimals=2)))
    # plt.savefig(local_path + "flanks_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def flank_plotting_reads(df, bct, appended):
    slope1, intercept1, _, _, _ = linregress(df['inside_array'], (df['right_flank'] + df['left_flank']))

    plt.title("Inside array count vs. flank counts for READ_LENGTH=" + str(read_length) + "\nSlope=" + str(np.round(slope1, decimals=2)) +
              ", intercept=" + str(np.round(intercept1, decimals=2)) + "\n" + bct)
    plt.ylabel("Flank count sum (right + left flanks)")
    plt.xlabel("Inside array count (inside_array)")
    plt.scatter(df['inside_array'], (df['right_flank'] + df['left_flank']))
    # plt.show()
    # plt.savefig(local_path + "flank_IA_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def flank_plotting_reads_driver():
    for (df, bct, appended) in [(df_n11, "Homozygous -1/-1", "n11"), (df_0n1, "Heterozygous 0/-1", "0n1"),
                                (df_00, "Homozygous 0/0", "00"), (df_01, "Heterozygous 0/1", "01"),
                                (df_11, "Homozygous 1/1", "11")]:
        flank_plotting_reads(df, bct, appended)

def inside_array_figures():
    for (df, bct, appended) in [(df_n11, "Homozygous -1/-1", "n11"), (df_0n1, "Heterozygous 0/-1", "0n1"),
                                (df_00, "Homozygous 0/0", "00"), (df_01, "Heterozygous 0/1", "01"),
                                (df_11, "Homozygous 1/1", "11")]:
        inside_array_plotting(df, bct, appended, np.mean(df['inside_array']), np.median(df['inside_array']))

def inside_array_plotting(df, bct, appended, df_mean, df_median):
    plt.hist((df['inside_array']), bins=20)
    plt.xlabel("Inside array count")
    plt.ylabel("Frequency")
    # plt.title("Histogram plot of inside_array count for READ_LENGTH=" + str(read_length) + "\n" + bct + "\nMean=" + str(np.round(df_mean, decimals=2)) + " | median=" + str(np.round(df_median, decimals=2)))
    # plt.savefig(local_path + "inside_array_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def pattern_ratio_plot(df, bct, appended):
    plt.title("Histogram of (array_length / pattern_size) READ_LENGTH=" + str(read_length) + "\n" + bct)
    plt.ylabel("Frequency")
    plt.xlabel("Ratio: (array_length / pattern_size)")
    # plt.hist(df['array_length'] / df['pattern_size'], bins=20)
    plt.hist(df['copy_number'], bins=20)
    plt.show()
    plt.clf()
    # plt.savefig(local_path + "flank_IA_" + appended + "cne" + str(copy_number_exclude) + "_" + str(read_length) + ".png")
    # plt.clf()

def pattern_ratio_plot_driver():
    for (df, bct, appended) in [(df_n11, "Homozygous -1/-1", "n11"), (df_0n1, "Heterozygous 0/-1", "0n1"),
                                (df_00, "Homozygous 0/0", "00"), (df_01, "Heterozygous 0/1", "01"),
                                (df_11, "Homozygous 1/1", "11")]:
        pattern_ratio_plot(df, bct, appended)

# inside_array_figures()
# flank_plotting_reads_driver()
# flank_figures()
# mean_and_std_figures()
# std_df.csv

# pattern_ratio_plot_driver()

# plt.scatter(df_00['inside_array'] / df_00['array_length'], np.log(df_00['ratio']))
def f(x, a, b, c):
    return a * (x ** 2) + b * x + c

from scipy.optimize import curve_fit

# df_00 = df_n11 #

popt, _ = curve_fit(f, np.log(df_00['pattern_size'] * df_00['ratio'] + 1), df_00['array_length'])

x_data = np.log(df_00['pattern_size'] * df_00['ratio'] + 1)
x_data2 = np.array(np.arange(4, 12, 0.1))
plt.scatter(x_data, df_00['array_length'])
plt.title("Prediction of array length by ratio and pattern size\nRead length = " + str(read_length))
plt.xlabel("Predictor = log(pattern_size * ratio)")
plt.ylabel("Prediction = array length")
print("Curve fitted parameters", popt)
plt.plot(x_data2,
         f(x_data2, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.show()

exit(0)

plt.hist(df_00['copy_number'], bins=20)
plt.title("0/0 " + str(np.array(df_00['copy_number']).mean()) +
            "\n0/1 " + str(np.array(df_01['copy_number']).mean()) +
            "\n1/1 " + str(np.array(df_11['copy_number']).mean()) +
            "\n0/-1 " + str(np.array(df_0n1['copy_number']).mean()) +
            "\n-1/-1 " + str(np.array(df_n11['copy_number']).mean()))