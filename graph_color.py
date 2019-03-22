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

def pandas_analyze(path, genotype_filter="", printE=0, title="", force_exclude=0, force_stop=0):
    dataframe = pd.read_csv(path)
    if force_stop != 0:
        dataframe = dataframe.iloc[:force_stop]
    dataframe = dataframe.loc[dataframe['filter'] == 'S']
    if len(genotype_filter) > 0:
        dataframe = dataframe.loc[dataframe['genotype'] == genotype_filter]
    dataframe = dataframe.loc[dataframe['right_flank'] > 0]
    dataframe = dataframe.loc[dataframe['left_flank'] > 0]
    dataframe['ratio'] = dataframe['inside_array'] / (dataframe['right_flank'] + dataframe['left_flank'])

    if force_exclude != 0:
        dataframe = dataframe.loc[dataframe['array_length'] < force_exclude]
    else:
        dataframe = dataframe.loc[dataframe['array_length'] < 10000]

    # inside_array = np.array(dataframe.iloc[:]["inside_array"])
    # right_flank = np.array(dataframe.iloc[:]["right_flank"])
    # left_flank = np.array(dataframe.iloc[:]["left_flank"])
    # array_length = np.array(dataframe.iloc[:]["array_length"])

    scaled_patterns = dataframe['pattern_size']
    scaled_patterns = scale_to_range(scaled_patterns, 0, 39 - 15)

    return dataframe, scaled_patterns

csv_filename = "sim3_100bp.csv"
df_00, colors1 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/0', force_exclude=5000, force_stop=1000)
df_01, colors2 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/1', force_exclude=5000, force_stop=1000)
df_11, colors3 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='43466', force_exclude=5000, force_stop=1000)
df_n01, colors4 = pandas_analyze("/Users/mkorovkin/Desktop/marzd/" + csv_filename,
                                genotype_filter='0/-1', force_exclude=5000, force_stop=1000)



trace1 = go.Scatter(
    x=df_00['array_length'],
    y=df_00['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=colors1,
        colorscale='Greens',
        showscale=True
    )
)
trace2 = go.Scatter(
    x=df_01['array_length'],
    y=df_01['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=colors2,
        colorscale='Blues',
        showscale=True
    )
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
    )
)
trace4 = go.Scatter(
    x=df_n01['array_length'],
    y=df_n01['ratio'],
    mode='markers',
    marker=dict(
        size=16,
        color=colors4,
        colorscale='Greys',
        showscale=True
    )
)
data = [trace1, trace2, trace3, trace4]

py.plot(data, filename='scatter-plot-with-colorscale')