import pandas as pd
import numpy as np
from scipy.stats import linregress
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt

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
    dataframe['ratio'] = np.array(dataframe['inside_array']) / 100.0 # / (dataframe['right_flank'] + dataframe['left_flank'])

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

def pandas_analyze_new(path, READ_LENGTH, flank_map, genotype_filter="", printE=0, title="", force_exclude=0,
                       force_stop=0, drop_lowest_percent_pattern_length=0.0, drop_lowest_percent_array_length=0.0,
                       copy_number_exclusion=0, pattern_array_length_ratio_exclude=0, force_flank=True):
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

    if pattern_array_length_ratio_exclude > 0:
        dataframe = dataframe.loc[dataframe['pattern_size'] > 0]
        dataframe['pattern_al'] = dataframe['array_length'] / dataframe['pattern_size']
        dataframe['pattern_al'] = dataframe.loc[dataframe['pattern_al'] < pattern_array_length_ratio_exclude]

    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]

    if copy_number_exclusion > 0:
        dataframe = dataframe.loc[dataframe['copy_number'] < copy_number_exclusion]

    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]
    if not force_flank:
        dataframe['ratio'] = dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank'])
    else:
        dataframe['ratio'] = dataframe['inside_array'] / flank_map[READ_LENGTH]#(dataframe['right_flank'] + dataframe['left_flank'])

    scaled_patterns = dataframe['pattern_size']
    scaled_patterns = scale_to_range(scaled_patterns, 0, 39 - 15)

    s, i, _, _, stderr = linregress(dataframe['array_length'], dataframe['ratio'])
    Alist = np.arange(READ_LENGTH * 3, READ_LENGTH * 50, READ_LENGTH)
    regression_array = np.array([(s * x + i) for x in Alist])

    error_bars = error_bars_around_dataset(dataframe, READ_LENGTH * 3, READ_LENGTH * 50, READ_LENGTH)

    return dataframe, np.log(5 + scaled_patterns), error_bars, regression_array

def get_mean_std(file, Alist, read_length, ins):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return np.array([(s * A + i) for A in ins]), np.array([(s1 * A + i1) for A in ins])

def show_hg_file2(drop_lowest_percent_pattern_length, flank_value, force_flank=True):
    file_path = "/Users/mkorovkin/Desktop/marzd/HG007counts.csv"

    df = pd.read_csv(file_path)

    df = df.loc[df['left_flank'] > 0]
    df = df.loc[df['right_flank'] > 0]
    df = df.loc[df['array_length'] > 0]
    df = df.loc[df['inside_array'] > 0]
    df = df.loc[df['array_length'] < 5000]
    df = df.loc[df['TRID'].str.contains("S")]

    if not force_flank:
        df['ratio'] = df['inside_array'] / (df['left_flank'] + df['right_flank'])
    else:
        df['ratio'] = df['inside_array'] / flank_value # (df['left_flank'] + df['right_flank'])

    df = df.loc[df['ratio'] < 200]

    file_path2 = "/Users/mkorovkin/Desktop/marzd/refset_full.csv"

    df2 = pd.read_csv(file_path2)
    df2.drop(['chrom', 'start', 'end'], axis=1)
    df2 = df2.loc[df2['indistinguishable'] == "S"]
    df2['TRID'] = df2['TRID'].apply(lambda x: str(x) + "_S")

    df = pd.merge(df, df2, on='TRID')

    x = df['array_length']

    if drop_lowest_percent_pattern_length > 0 and drop_lowest_percent_pattern_length < 1:
        df = upper_percentile_ps(df, drop_lowest_percent_pattern_length)
        df = df.reset_index()

    scaled_patterns = df['pattern_size']
    scaled_patterns = np.log(scale_to_range(scaled_patterns, 0, 39 - 15) + 1)

    #pat_list = list()

    #for pat in range(len(df['pattern_size'])):
    #    t = (df.loc[pat]['pattern_size'], df.loc[pat]['copy_number'])
    #    pat_list.append(t)

    return x, df['ratio'], scaled_patterns, df['pattern_size']

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

def run_script(read_length, smooth_std, excludeperc_ps=-0.7, excludeperc_arraylength=-0.9, copy_number_exclude=8, pattern_array_length_ratio_exclude2=10):
    # read_length = 148
    Alist2 = np.arange(read_length * 3, read_length * (40 - 7), read_length)
    Alist = np.arange(read_length * 3, read_length * 50, read_length)
    coverage = np.floor((355 - 70) / 2)

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

    #excludeperc_ps = -0.7 # 0.4
    #excludeperc_arraylength = -0.9

    flank_map = {100: 102.0,
                 148: 103.0,
                 250: 194.0}

    xHG, ratioHG, colorsHG, patternsHG = show_hg_file2(excludeperc_ps, flank_map[read_length])

    alist_standard, stdlist_standard = score_df_for_std_from_file(read_length, 100)
    _, _, ypoints = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, alist_standard, force_exclude=maxax, static=(True, read_length))

    #copy_number_exclude = 8
    #pattern_array_length_ratio_exclude2 = 10

    df_00, color00, arraylength_data00, regression_arr00 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length, flank_map,
                                                                              genotype_filter='0/0', force_exclude=maxax,
                                                                              drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                              drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                              pattern_array_length_ratio_exclude=pattern_array_length_ratio_exclude2,
                                                                              copy_number_exclusion=copy_number_exclude
                                                                              )
    df_0n1, color0n1, arraylength_data0n1, regression_arr0n1 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length, flank_map,
                                                                                  genotype_filter='0/-1', force_exclude=maxax,
                                                                                  drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                                  drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                                  pattern_array_length_ratio_exclude=pattern_array_length_ratio_exclude2,
                                                                                    copy_number_exclusion=copy_number_exclude)
    df_01, color01, arraylength_data01, regression_arr01 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length, flank_map,
                                                                              genotype_filter='0/1', force_exclude=maxax,
                                                                              drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                              drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                              pattern_array_length_ratio_exclude=pattern_array_length_ratio_exclude2,
                                                                              copy_number_exclusion=copy_number_exclude)
    df_11, color11, arraylength_data11, regression_arr11 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length, flank_map,
                                                                              genotype_filter='43466', force_exclude=maxax,
                                                                              drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                              drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                              pattern_array_length_ratio_exclude=pattern_array_length_ratio_exclude2,
                                                                              copy_number_exclusion=copy_number_exclude)
    df_n11, colorn11, arraylength_datan11, regression_arrn11 = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length, flank_map,
                                                                                  genotype_filter='1', force_exclude=maxax,
                                                                                  drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                                  drop_lowest_percent_array_length=excludeperc_arraylength,
                                                                                  pattern_array_length_ratio_exclude=pattern_array_length_ratio_exclude2,
                                                                                    copy_number_exclusion=copy_number_exclude)

    fitted_means_expected, fitted_stds = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length, Alist2)

    show00 = True
    show01 = True
    show11 = True
    show0n1 = True
    show_general00 = (True, True)
    show_expected = (True, False)
    show_whole_dataset_best_fit = (False, False)

    calpha = 1

    def get_trace_0():
        trace0_0 = go.Scatter(
            x=Alist2,
            y=fitted_means_expected,
            mode='lines',
            marker=dict(
                size=16,
                color=[0, 0],
                colorscale='Greys',
            ),
            name="(1) Expected regression line based on simulations"
        )

        trace0_1 = go.Scatter(
            x=Alist2,
            y=fitted_means_expected + fitted_stds,
            mode='lines',
            marker=dict(
                size=16,
                color=[0, 0],
                colorscale='Greys',
            ),
            name="(1)+Std based on simulations"
        )

        trace0_2 = go.Scatter(
            x=Alist2,
            y=fitted_means_expected - fitted_stds,
            mode='lines',
            marker=dict(
                size=16,
                color=[0, 0],
                colorscale='Greys',
            ),
            name="(1)-Std based on simulations"
        )

        trace0 = go.Scatter(
                x=xHG,
                y=ratioHG,
                mode='markers',
                marker=dict(
                    size=16,
                    color=colorsHG,
                    colorscale='Reds',
                    showscale=True
                ),
                name="HG007 Points",
                hovertext=patternsHG
            )

        return trace0, trace0_0, trace0_1, trace0_2

    trace0, trace0_0, trace0_1, trace0_2 = get_trace_0()
    mAlist = np.arange(read_length * 3, read_length * 50, read_length)

    ########### TRACE of 0/0
    trace1 = go.Scatter(
        x=df_00['array_length'],
        y=df_00['ratio'],
        mode='markers',
        marker=dict(
            size=16,
            color=color00,
            colorscale='Reds',
            showscale=True
        ),
        name="Homozygous (0/0) simulated points",
        hovertext=df_00['pattern_size']
    )
    trace1_copynumber = go.Scatter(
        x=df_00['array_length'],
        y=df_00['ratio'],
        mode='markers',
        marker=dict(
            size=16,
            color=color00,
            colorscale='Reds',
            showscale=True
        ),
        name="Homozygous (0/0) simulated points (copy number)",
        hovertext=df_00['copy_number']
    )
    trace1_0 = go.Scatter(
        x=mAlist,
        y=regression_arr00,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(0/0) - Expected regression line"
    )

    trace1_1 = go.Scatter(
        x=mAlist,
        y=regression_arr00 + arraylength_data00,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(0/0) - Std. Lines (+) reg. line"
    )

    trace1_2 = go.Scatter(
        x=mAlist,
        y=regression_arr00 - arraylength_data00,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(0/0) - Std. Lines (-) reg. line"
    )

    ########### TRACE of 0/-1
    trace2 = go.Scatter(
                x=df_0n1['array_length'],
                y=df_0n1['ratio'],
                mode='markers',
                marker=dict(
                    size=16,
                    color=color0n1,
                    colorscale='Greens',
                    showscale=True
                ),
                name="Heterozygous (0/-1) simulated points",
                hovertext=df_0n1['pattern_size']
            )
    trace2_0 = go.Scatter(
        x=mAlist,
        y=regression_arr0n1,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(0/-1) - Expected regression line"
    )

    ########### TRACE of 0/1
    trace3 = go.Scatter(
                x=df_01['array_length'],
                y=df_01['ratio'],
                mode='markers',
                marker=dict(
                    size=16,
                    color=color01,
                    colorscale='Blues',
                    showscale=True
                ),
                name="Heterozygous (0/1) simulated points",
                hovertext=df_01['pattern_size']
            )
    trace3_0 = go.Scatter(
        x=mAlist,
        y=regression_arr01,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(0/1) - Expected regression line"
    )

    ########### TRACE of 1/1
    trace4 = go.Scatter(
                x=df_11['array_length'],
                y=df_11['ratio'],
                mode='markers',
                marker=dict(
                    size=16,
                    color=color11,
                    colorscale='Electric',
                    showscale=True
                ),
                name="Homozygous (1/1) simulated points",
                hovertext=df_11['pattern_size']
            )
    trace4_0 = go.Scatter(
        x=mAlist,
        y=regression_arr11,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(1/1) - Expected regression line"
    )

    ########### TRACE of -1/-1
    trace5 = go.Scatter(
                x=df_n11['array_length'],
                y=df_n11['ratio'],
                mode='markers',
                marker=dict(
                    size=16,
                    color=colorn11,
                    colorscale='Earth',
                    showscale=True
                ),
                name="Homozygous (-1/-1) simulated points",
                hovertext=df_n11['pattern_size']
            )
    trace5_0 = go.Scatter(
        x=mAlist,
        y=regression_arrn11,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(-1/-1) - Expected regression line"
    )

    def get_written_std(read_length, smooth=False):
        input_path = ""
        if read_length == 100:
            input_path = "/Users/mkorovkin/Desktop/graphs_marzd/std_df_cne8_100.csv"
        elif read_length == 148:
            input_path = "/Users/mkorovkin/Desktop/graphs_marzd/std_df_cne8_148.csv"
        elif read_length == 250:
            input_path = "/Users/mkorovkin/Desktop/graphs_marzd/std_df_cne8_250.csv"
        temp_df = pd.read_csv(input_path)

        stda = temp_df['stda'] * 1.96
        stdx = temp_df['stdx']

        if not smooth:
            return stda, stdx
        else:
            slope1, intercept1, _, _, _ = linregress(stdx, stda)
            return [(x * slope1 + intercept1) for x in Alist], Alist

    stda, stdx = get_written_std(read_length, smooth_std)
    fitmeans, _ = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length, stdx)

    trace_stdp = go.Scatter(
        x=stdx,
        y=fitmeans + stda,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(2) (+ Grouped std) - Expected regression line"
    )

    trace_std = go.Scatter(
        x=stdx,
        y=fitmeans * (3.8 - 0.5) / 3.8,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(2) - Expected regression line"
    )

    trace_stdm = go.Scatter(
        x=stdx,
        y=fitmeans - stda,
        mode='lines',
        marker=dict(
            size=16,
            colorscale='Greys',
            showscale=True
        ),
        name="(2) (- Grouped std) - Expected regression line"
    )

    layout = go.Layout(title="READ_LENGTH=" + str(read_length) + " | coverage=" + str(coverage) + " | array lengths from [" + str(read_length * 3) + ":" + str(read_length * 50) + "]",
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

    data = [trace0, trace1, trace2, trace3, trace4, trace5, #trace1_copynumber,
            # trace1_1, trace1_2,
            trace1_0, trace2_0, trace3_0, trace4_0, trace5_0,
            # trace0_0, trace0_1, trace0_2,
            trace_std, trace_stdm, trace_stdp]#, trace1, trace1_0, trace1_1, trace1_2]

    py.plot(go.Figure(data, layout), filename='scatter-plot-with-colorscale')

# Clarify: key things to focus on
    # When to go to genome; what kind of model to use with the genome