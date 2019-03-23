import pandas as pd
import numpy as np
from scipy.stats import linregress
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

def scale_to_range(list_, min_, max_):
    list_ = np.array(list_, dtype=np.float32)
    min_sp = list_.min()
    max_sp = list_.max()
    list_ = list_ - min_sp
    list_ = list_ / max_sp
    list2 = (1 - list_) * min_ / max_
    list_ = (list_ + list2) * max_
    return np.array(np.round(list_, decimals=0), dtype=np.int32)

# print(scale_to_range([1, 2, 3, 4, 5, 6, 7, 8], 0, 39-15))

def invert_colors(array):
    minx = array.min()
    maxx = array.max()

    return array * -1 + maxx + minx

def pandas_analyze(path, genotype_filter="", printE=0, title="", force_exclude=0, force_stop=0):
    dataframe = pd.read_csv(path)
    if force_stop != 0:
        dataframe = dataframe.iloc[:force_stop]

    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    else:
        dataframe = dataframe.loc[dataframe['array_length'] < 15000]

    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]
    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]
    dataframe['ratio'] = dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank'])

    # inside_array = np.array(dataframe.iloc[:]["inside_array"])
    # right_flank = np.array(dataframe.iloc[:]["right_flank"])
    # left_flank = np.array(dataframe.iloc[:]["left_flank"])
    # array_length = np.array(dataframe.iloc[:]["array_length"])

    scaled_patterns = dataframe['pattern_size']
    scaled_patterns = scale_to_range(scaled_patterns, 0, 39 - 15)

    return dataframe, scaled_patterns

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
            return [minx, maxx], [0.00873 * minx + -1.57248, 0.00873 * maxx + -1.57248], np.array([(0.00873 * x + -1.57248) for x in xs])
        elif static[1] == 148:
            return [minx, maxx], [0.00669 * minx + -2.03295, 0.00669 * maxx + -2.03295], np.array([(0.00669 * x + -2.03295) for x in xs])
        elif static[1] == 250:
            return [minx, maxx], [0.00491 * minx + -2.68229, 0.00491 * maxx + -2.68229], np.array([(0.00491 * x + -2.68229) for x in xs])
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

def score_df_for_std(bp, coverage):
    eq = (0, 0)
    if bp == 100:
        READ_LENGTH = 100
        Alist = np.arange(READ_LENGTH * 3, READ_LENGTH * 50, READ_LENGTH)
        std_list = list()
        if coverage == 70:
            std_list = np.array(
            [0.529, 1.0, 1.45, 1.49, 2.136, 2.136, 2.136, 3.038, 4.038, 4.077, 4.998, 5.272, 6.292, 5.841, 7.468, 7.84,
             9.408, 8.506, 9.114, 9.604, 10.016, 10.937, 9.114, 11.584, 13.21, 13.093, 12.838, 15.072, 15.621, 17.013,
             17.248, 16.836, 15.504, 25.068, 17.718, 19.6, 17.954, 22.697, 16.66, 22.226, 21.991, 22.658, 22.011,
             22.089, 25.519, 23.559, 25.225])
        elif coverage == 355:
            std_list = np.array([0.235, 0.431, 0.686, 0.706, 0.764, 0.921, 1.137, 1.196, 1.411, 1.901, 2.038, 1.94, 1.96, 3.175, 2.666,
                  3.018, 3.41, 3.273, 3.41, 3.94, 4.175, 4.528, 5.037, 4.371, 5.547, 5.841, 5.41, 6.35, 5.684, 6.037,
                  6.488, 7.036, 6.958, 6.88, 6.938, 7.526, 7.781, 7.683, 8.624, 9.016, 8.271, 9.643, 9.741, 8.173,
                  8.996, 10.937, 8.82])
        else:
            std_list = np.array(
            [0.382, 0.715, 1.068, 1.098, 1.45, 1.529, 1.637, 2.117, 2.724, 2.989, 3.518, 3.606, 4.126, 4.508, 5.067,
             5.429, 6.409, 5.89, 6.262, 6.772, 7.095, 7.732, 7.076, 7.977, 9.379, 9.467, 9.124, 10.711, 10.653, 11.525,
             11.868, 11.936, 11.231, 15.974, 12.328, 13.563, 12.867, 15.19, 12.642, 15.621, 15.131, 16.15, 15.876,
             15.131, 17.258, 17.248, 17.023])
        return Alist, std_list#[minx, maxx], np.array([minx * eq[0] + eq[1], maxx * eq[0] + eq[1]])
    if bp == 148:
        READ_LENGTH = 148
        Alist = np.arange(READ_LENGTH * 3, READ_LENGTH * 50, READ_LENGTH)
        std_list = list()
        if coverage == 70:
            std_list = np.array([0.627, 1.058, 1.156, 1.372, 1.744, 2.47, 3.156, 3.822, 3.842, 4.743, 5.547, 6.056, 7.37, 6.311, 7.82,
                  8.624, 7.134, 9.996, 8.742, 10.408, 11.584, 10.466, 12.818, 16.327, 14.406, 13.367, 14.661, 13.544,
                  15.954, 14.23, 14.367, 19.071, 16.915, 19.737, 15.876, 18.463, 26.656, 21.364, 19.267, 21.658, 21.501,
                  20.325, 24.304, 23.089, 23.716, 29.733, 25.598])
        elif coverage == 355:
            std_list = np.array(
            [0.235, 0.372, 0.49, 0.608, 0.706, 1.0, 1.274, 1.509, 1.705, 1.882, 2.058, 2.47, 2.509, 2.822, 2.999, 2.842,
             3.293, 3.9, 3.646, 3.802, 4.371, 4.978, 4.978, 5.41, 4.978, 5.468, 5.351, 5.9, 5.723, 6.488, 7.311, 7.82,
             7.938, 7.252, 7.624, 7.84, 7.938, 8.271, 8.095, 7.33, 9.232, 8.722, 9.506, 11.211, 10.604, 10.27, 9.114])
        else:
            std_list = np.array([0.431, 0.715, 0.823, 0.99, 1.225, 1.735, 2.215, 2.666, 2.773, 3.312, 3.802, 4.263, 4.939, 4.567, 5.41,
                  5.733, 5.214, 6.948, 6.194, 7.105, 7.977, 7.722, 8.898, 10.868, 9.692, 9.418, 10.006, 9.722, 10.839,
                  10.359, 10.839, 13.446, 12.426, 13.495, 11.75, 13.152, 17.297, 14.818, 13.681, 14.494, 15.366, 14.524,
                  16.905, 17.15, 17.16, 20.002, 17.356])
        return Alist, std_list
    if bp == 250:
        READ_LENGTH = 250
        Alist = np.arange(READ_LENGTH * 3, READ_LENGTH * 50, READ_LENGTH)
        std_list = list()
        if coverage == 70:
            std_list = np.array(
        [0.451, 0.706, 1.313, 1.529, 2.136, 2.901, 3.41, 4.096, 3.45, 5.116, 5.233, 5.802, 5.978, 8.036, 8.252, 8.879,
         9.055, 12.76, 8.271, 9.232, 11.407, 12.093, 12.27, 11.113, 12.407, 14.837, 15.641, 13.348, 16.092, 16.856,
         17.758, 21.521, 19.482, 19.502, 19.835, 21.364, 17.385, 19.286, 19.796, 21.815, 20.031, 21.815, 23.167, 24.01,
         22.638, 25.97, 26.225])
        elif coverage == 355:
            std_list = np.array(
        [0.216, 0.333, 0.51, 0.706, 0.784, 1.058, 1.274, 1.646, 1.803, 2.332, 2.097, 2.43, 2.881, 2.999, 3.43, 3.802,
         3.273, 4.234, 4.038, 4.273, 4.724, 4.469, 5.233, 5.86, 5.214, 5.939, 5.704, 5.86, 7.036, 6.84, 5.645, 6.566,
         7.428, 8.252, 7.213, 8.114, 8.683, 7.115, 8.898, 7.879, 8.428, 9.408, 10.251, 9.898, 10.78, 10.447, 11.348 ])
        else:
            std_list = np.array(
        [0.333, 0.519, 0.911, 1.117, 1.46, 1.98, 2.342, 2.871, 2.626, 3.724, 3.665, 4.116, 4.43, 5.517, 5.841, 6.341,
         6.164, 8.497, 6.154, 6.752, 8.065, 8.281, 8.751, 8.487, 8.81, 10.388, 10.672, 9.604, 11.564, 11.848, 11.701,
         14.043, 13.455, 13.877, 13.524, 14.739, 13.034, 13.201, 14.347, 14.847, 14.23, 15.611, 16.709, 16.954, 16.709,
         18.208, 18.787])
        return Alist, std_list
    return
    #(0.00531795744680851, -1.5549748226950353) <- for C=70
    #(0.0022562765957446804, -0.8927053191489343) <- for C=355
    #(0.0037871170212765964, -1.2238400709219857) <- for MEAN

read_length = 100
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

df_00, colors1 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/0', force_exclude=maxax)#, filter_SI=True)
df_01, colors2 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/1', force_exclude=maxax)#, filter_SI=True)
df_11, colors3 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='43466', force_exclude=maxax)#, filter_SI=True)
df_n01, colors4 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/-1', force_exclude=maxax)#, filter_SI=True)

xsin, stds = score_df_for_std(read_length, coverage)
_, _, ypoints = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, xsin, force_exclude=maxax)
_, _, ypoints2 = score_df_nofilter("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, xsin, force_exclude=maxax)
_, _, ypoints3 = score_df("/Users/mkorovkin/Desktop/marzd/" + csv_filename, minax, maxax, xsin, force_exclude=maxax, static=(True, read_length))
calpha = 1

trace0 = go.Scatter(
    x=xsin,
    y=ypoints,
    mode='lines',
    marker=dict(
        size=16,
        color=[0, 0],
        colorscale='Greys',
    ),
    name="(1) Line of best fit (based on 0/0 genotype)"
)
trace00 = go.Scatter(
    x=xsin,
    y=ypoints2,
    mode='lines',
    marker=dict(
        size=16,
        color=[0, 0],
        colorscale='Greys',
    ),
    name="(2) Line of best fit of the whole dataset"
)
trace01 = go.Scatter(
    x=xsin,
    y=ypoints3,
    mode='lines',
    marker=dict(
        size=16,
        color=[0, 0],
        colorscale='Greys',
    ),
    name="(3) Expected regression line"
)
trace0_1 = go.Scatter(
    x=xsin,
    y=(ypoints - stds),
    mode='lines',
    marker=dict(
        size=16,
        color=[0, 0],
        colorscale='Reds',
    ),
    name="* Lower standard deviation line"
)
trace0_2 = go.Scatter(
    x=xsin,
    y=(ypoints + stds),
    mode='lines',
    marker=dict(
        size=16,
        color=[0, 0],
        colorscale='Reds',
    ),
    name="* Upper standard deviation line"
)
trace1 = go.Scatter(
    x=df_00['array_length'],
    y=df_00['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=invert_colors(colors1),
        colorscale='Greens',
        showscale=True
    ),
    name="Homozygous: 0/0",
    hovertext=df_00['pattern_size']
)
trace2 = go.Scatter(
    x=df_01['array_length'],
    y=df_01['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=invert_colors(colors2),
        colorscale='Blues',
        showscale=True
    ),
    name="Heterozygous: 0/1",
    hovertext=df_01['pattern_size']
)
trace3 = go.Scatter(
    x=df_11['array_length'],
    y=df_11['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=colors3,
        colorscale='Reds',
        showscale=True
    ),
    name="Homozygous: 1/1",
    hovertext=df_11['pattern_size']
)
trace4 = go.Scatter(
    x=df_n01['array_length'],
    y=df_n01['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=colors4,
        colorscale='Viridis',
        showscale=True
    ),
    name="Heterozygous: 0/-1",
    hovertext=df_n01['pattern_size']
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
data = [trace4, trace1, trace2, trace3, trace0, trace00, trace01, trace0_1, trace0_2]

py.plot(go.Figure(data, layout), filename='scatter-plot-with-colorscale')

# KEY NOTES: higher pattern lengths tend to cause more deviation from the trend; depends on whether it's a gain or a loss
# 5% outliers around the line
# read up on the statistics/standard deviation
# histogram of where the points are (density fo points)
# in bin sizes: piecewise linear regression
    # means in theory
# count fragment ratios; confirm the reasoning behind association of fragment length
# remove buffers on the inside (25% buffer)
