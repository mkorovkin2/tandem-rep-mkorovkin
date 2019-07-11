import pandas as pd
import numpy as np

df = pd.read_csv("/Users/mkorovkin/Desktop/marzd/HG007counts.csv")
df = df.loc[df['right_flank'] > 0]
df = df.loc[df['left_flank'] > 0]
df = df.loc[df['inside_array'] > 0]
df['ratio'] = df['inside_array'] / (df['right_flank'] + df['left_flank'])
print(df)
df.to_csv("/Users/mkorovkin/Desktop/marzd/HG007counts_mfhy.csv")