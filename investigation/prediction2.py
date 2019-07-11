import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

def label_row(row):
   if row['genotype'] == '1':
      return -1
   elif row['genotype'] == '43466':
      return 1
   elif row['genotype'] == '0/1':
       return 0.5
   elif row['genotype'] == '0/-1':
       return -0.5
   elif row['genotype'] == '0/0':
       return 0
   else:
       return 0

def label_correct_negative(row, fit_mean, fit_std):
    y = row['ratio']
    arr = row['array_length']
    pat = row['pattern_size']

    x = arr - pat

    expected_y = fit_mean[0] * x + fit_mean[1]
    expected_std = fit_std[0] * x + fit_std[1]
    expected_range = (expected_y - expected_std, expected_y + expected_std)

    if y > expected_range[1]:
        return 1
    elif y < expected_range[0]:
        return -1
    else:
        return 0

def label_indistinguish_negative(row, fit_mean, fit_std):
    y = row['ratio']
    x = row['array_length']

    expected_y = fit_mean[0] * x + fit_mean[1]
    expected_std = fit_std[0] * x + fit_std[1]
    expected_range = (expected_y - expected_std, expected_y + expected_std)

    if y < expected_range[0] or y > expected_range[1]:
        return 1
    else:
        return 0

def get_stats(file, Alist, read_length):
    df = pd.read_csv(file)[:30]
    means = np.array(df["mA" + str(read_length)])
    stds = np.array(df["sA" + str(read_length)])

    s, i, _, _, _ = linregress(Alist, means)
    s1, i1, _, _, _ = linregress(Alist, stds)

    return (s, i), (s1, i1)

def percent_error(row, mean_inside_array, abs_or_not):
    x0 = row['inside_array']

    if abs_or_not:
        return abs(1 - (x0 / mean_inside_array))
    else:
        return 1 - (x0 / mean_inside_array)

# Calculate necessary difference, based on standard deviation bars, for a given array length
def calculate_statistical_difference(x, mean, std):
    beta1 = (mean[0] + std[0])
    beta0 = mean[1] + std[1]
    gamma1 = (mean[0] - std[0])
    gamma0 = mean[1] + std[1]

    return ((beta1 - gamma1) * x + (beta0 - gamma0)) / gamma1 # Difference in base pairs

read_length = 148

csv_filename = ""
if read_length == 100:
    csv_filename = "sim3_100bp.csv"
elif read_length == 148:
    csv_filename = "sim2_148bp.csv"
elif read_length == 250:
    csv_filename = "sim1_250bp.csv"
else:
    exit(0)

df = pd.read_csv("/Users/mkorovkin/Desktop/marzd/" + csv_filename)

# Clean the data
df = df.loc[df['filter'] == "S"]
df = df.loc[df['right_flank'] > 0]
df = df.loc[df['left_flank'] > 0]
df = df.loc[df['inside_array'] > 0]

# Build new columns
df['label'] = df.apply(lambda row: label_row(row), axis=1)
df['ratio'] = np.array(df['inside_array']) / np.array(df['right_flank'] + df['left_flank'])

# Perform analysis of -1/-1
Alist2 = np.arange(read_length * 3, read_length * (40 - 7), read_length)
fit_mean, fit_std = get_stats("/Users/mkorovkin/Desktop/marzd/output40_new.csv", Alist2, read_length)

dfn = df.loc[df['label'] == -1]

dfn['predictable'] = dfn.apply(lambda row: label_correct_negative(row, fit_mean, fit_std), axis=1)
dfn['unpredictable'] = dfn.apply(lambda row: label_indistinguish_negative(row, fit_mean, fit_std), axis=1)

# Perform basic analysis on requirements for distinguishibility
Alist = np.arange(read_length * 3, read_length * 50, read_length)
basic_reg_line = np.array([(fit_mean[0] * x + fit_mean[1]) for x in Alist])
basic_reg_line_with_lower_std = np.array([(fit_mean[0] * x + fit_mean[1] - (fit_std[0] * x + fit_std[1])) for x in Alist])

line_difference = basic_reg_line - basic_reg_line_with_lower_std

percentage_difference = line_difference / basic_reg_line

plt.plot(Alist, percentage_difference)
plt.title("Array Length & Prediction of Statistical Error")
plt.ylabel("Proportion error between two ratios required for 95% Z-test")
plt.xlabel("Array length")
plt.show()

lbound = 500
ubound = 2000#lbound + 500
df2 = df.loc[df['array_length'] > lbound]
df2 = df2.loc[df['array_length'] < ubound]

sum_flanks = np.array(df2['right_flank']) + np.array(df2['left_flank'])
expected_flank_sum = sum_flanks.mean()

print("Flank sum", expected_flank_sum)
print("Expected ratio mean", np.array(df2['inside_array']).mean() / expected_flank_sum)
print("Ratio mean", np.array(df2['ratio']).mean())
print("Mean array length", np.array(df2['array_length']).mean())
print("---")

mean_inside_array = np.array(df2['inside_array']).mean()

df2['error'] = df2.apply(lambda row: percent_error(row, mean_inside_array, False), axis=1)
df2['abs_error'] = df2.apply(lambda row: percent_error(row, mean_inside_array, True), axis=1)

plt.hist(df2['error'], bins=40)
plt.title("Frequency of Percent Error for Data Subset A=[" + str(lbound) + ", " + str(ubound) + "]")
plt.ylabel("Frequency of percent error")
plt.xlabel("Percent Error")
plt.show()

print("Mean error", np.array(df2['error']).mean())
print("95th error percentile", np.percentile(df2['abs_error'], 95))
print("---")

print("Prop. above 25% error", len(df2.loc[df2['abs_error'] > 0.25]) / len(df2))
print("*   ", len(df2.loc[df2['abs_error'] > 0.25]), "/", len(df2))
df11 = df2.loc[df2['genotype'] == '43466']
df01 = df2.loc[df2['genotype'] == '0/1']
df00 = df2.loc[df2['genotype'] == '0/0']
dfn01 = df2.loc[df2['genotype'] == '0/-1']
dfn11 = df2.loc[df2['genotype'] == '1']
print("Prop. of 1/1 above 25% error", len(df11.loc[df11['abs_error'] > 0.25]) / len(df11))
print("Prop. of 0/1 above 25% error", len(df01.loc[df01['abs_error'] > 0.25]) / len(df01))
print("Prop. of 0/0 above 25% error", len(df00.loc[df00['abs_error'] > 0.25]) / len(df00))
print("Prop. of 1/-1 above 25% error", len(dfn01.loc[dfn01['abs_error'] > 0.25]) / len(dfn01))
print("Prop. of -1/-1 above 25% error", len(dfn11.loc[dfn11['abs_error'] > 0.25]) / len(dfn11))

# Calculate lbound/ubound requirements approximations for significance
#   print(fit_std[0] * 1000 + fit_std[1] + fit_mean[0] * 1000 + fit_mean[1])
#   print((fit_mean[0] * 1000 + fit_mean[1]) * 1.255)

print("---")
print("Required base pair difference between array lengths")
limlist = list()
for A in Alist:
    limit = np.ceil(calculate_statistical_difference(A, fit_mean, fit_std))
    print("A=" + str(A) + ": " + str(np.int32(limit)) + " base pairs")
    limlist.append(limit)

plt.plot(Alist, limlist)
plt.title("Required array length for statistical difference between array lengths")
plt.ylabel("Array base pair difference required")
plt.xlabel("Array length")
plt.show()

# Distribution of flanks
# Making assumptions; cannot retain information (filter array length) -> indicate that things are done
#