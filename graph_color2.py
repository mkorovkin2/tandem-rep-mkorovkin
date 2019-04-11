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

def pandas_analyze_new(path, READ_LENGTH, genotype_filter="", printE=0, title="", force_exclude=0, force_stop=0, drop_lowest_percent_pattern_length=0.0, drop_lowest_percent_array_length=0.0):
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

    s, i, _, _, stderr = linregress(dataframe['array_length'], dataframe['ratio'])
    Alist = np.arange(read_length * 3, read_length * 50, read_length)
    yeb = np.array([(s * x + i) for x in Alist])
    stderr = stderr * len(dataframe['ratio'])
    return dataframe, scaled_patterns, np.array(dataframe['array_length']), yeb, yeb + stderr * 1.96, yeb - stderr * 1.96

def get_mean_std(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return np.array([(s * A + i) for A in Alist]), np.array([(s1 * A + i1) for A in Alist])

def show_hg_file():
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

    x = df['array_length']

    return x, df['ratio']

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

    scaled_patterns = df['pattern_size']
    scaled_patterns = np.log(scale_to_range(scaled_patterns, 0, 39 - 15) + 1)

    pat_list = list()
    for pat in range(len(df['pattern_size'])):
        t = (df.loc[pat]['pattern_size'], df.loc[pat]['copy_number'])
        pat_list.append(t)

    return x, df['ratio'], scaled_patterns, df['pattern_size']

read_length = 148
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

excludeperc_ps = -0.7
excludeperc_arraylength = -0.9

xHG, ratioHG, colorsHG, patternsHG = show_hg_file2()

alist_standard, stdlist_standard = score_df_for_std_from_file(read_length, 100)
_, _, ypoints = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, alist_standard, force_exclude=maxax, static=(True, read_length))

df_00, colors1, arraylength_data00, yeb, yebp, yebm = pandas_analyze_new("/Users/mkorovkin/Desktop/marzd/" + csv_filename, read_length,
                                                                   genotype_filter='0/0', force_exclude=maxax,
                                                                   drop_lowest_percent_pattern_length=excludeperc_ps,
                                                                   drop_lowest_percent_array_length=excludeperc_arraylength)

fitted_means_expected, fitted_stds = get_mean_std("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length)

show00 = True
show01 = True
show11 = True
show0n1 = True
show_general00 = (True, True)
show_expected = (True, False)
show_whole_dataset_best_fit = (False, False)

colors1 = np.log(colors1 + 5)

# _, _, ypoints = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, xsin, force_exclude=maxax)
calpha = 1

trace0 = go.Scatter(
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

trace2 = go.Scatter(
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

def dump():
    trace4 = go.Scatter(
        x=arraylength_data00,
        y=df_00['ratio'],
        mode='markers',
        marker=dict(
            size=16,
            color=colors1,
            colorscale='Viridis',
        ),
        name="(2) 0/0 Homozygous regressive fit line"
    )

    trace1 = go.Scatter(
        x=Alist,
        y=yeb,
        mode='lines',
        marker=dict(
            size=16,
            color=colors1,
            colorscale='Viridis',
        ),
        name="(2) Regressive fit for homozygous dataset"
    )

    trace3 = go.Scatter(
        x=alist_standard,
        y=ypoints,
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="(3) Expected regression line based on dataset (last week's simulations)"
    )

    trace3_1 = go.Scatter(
        x=alist_standard,
        y=ypoints - stdlist_standard,
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="(3)-Std on Expected regression line"
    )

    trace3_2 = go.Scatter(
        x=alist_standard,
        y=ypoints + stdlist_standard,
        mode='lines',
        marker=dict(
            size=16,
            color=[0, 0],
            colorscale='Greys',
        ),
        name="(3)+Std on Expected regression line"
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

data = [trace2, trace0, trace0_1, trace0_2]

py.plot(go.Figure(data, layout), filename='scatter-plot-with-colorscale')

# Clarify: key things to focus on
    # When to go to genome; what kind of model to use with the genome