import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from scipy.stats import linregress
import tensorflow as tf

def simulation_set():
    hg_set = pd.read_csv("/Users/mkorovkin/Desktop/marzd/sim2_148bp.csv")

    hg_set.drop("start", axis=1)
    hg_set.drop("end", axis=1)

    hg_set = hg_set.loc[hg_set["filter"] == "S"]
    hg_set = hg_set.loc[hg_set["array_length"] < 15000]
    hg_set = hg_set.loc[hg_set["pattern_size"] > 150]
    hg_set = hg_set.loc[hg_set["left_flank"] > 0]
    hg_set = hg_set.loc[hg_set["right_flank"] > 0]
    hg_set = hg_set.loc[hg_set["inside_array"] > 0]

    # Keep flanks static
    lfm = hg_set.left_flank.mean()
    rfm = hg_set.right_flank.mean()
    # print("flanks", lfm, rfm)
    hg_set["left_flank"] = hg_set.left_flank.apply(lambda x: lfm)
    hg_set["right_flank"] = hg_set.right_flank.apply(lambda x: rfm)

    # Calculate ratios
    hg_set["ratio"] = hg_set.inside_array / (hg_set.left_flank + hg_set.right_flank)
    hg_set["array_over_pattern"] = hg_set.array_length / hg_set.pattern_size
    hg_set = hg_set.loc[hg_set["array_over_pattern"] < 10]

    # prepare statistical data
    sm = pd.read_csv("/Users/mkorovkin/Desktop/output_statistics_00_new.csv")
    xs = np.arange(148 * 3, 148 * 40, 148)
    ms = sm.mA148

    slope1, intercept1, _, _, _ = linregress(xs, ms)

    # analyze total_set
    hg_set['inside_length_ratio'] = hg_set.inside_array / hg_set.array_length
    hg_set['error'] = (hg_set.array_length.apply(lambda x: slope1 * x + intercept1) - hg_set.ratio) / hg_set.ratio

    zin = {
        "-1/-1": "43466",
        "1/1": "1",
        "0/1": "0/1",
        "0/-1": "0/-1",
        "0/0": "0/0",
    }
    nex = {
        "43466": -2,
        "1": 2,
        "0/1": 1,
        "0/-1": -1,
        "0/0": 0,
    }

    hg_set['label'] = hg_set.genotype.apply(lambda x: nex[x])

    return hg_set

def linear_regression():
    x1 = tf.placeholder(tf.float32, name="x1")
    x2 = tf.placeholder(tf.float32, name="x2")
    y = tf.placeholder(tf.float32, name="y")

    with tf.variable_scope('lreg') as scope:
        b1 = tf.Variable(np.random.normal(), name="b1")# / 10 + 0.0001, name="b1")# 0.0001, name="b1")
        b2 = tf.Variable(np.random.normal(), name="b2")# / 10 + 0.0001, name="b2")# 0.0001, name="b2")
        b0 = tf.Variable(np.random.normal(), name="b0")# / 10 + 0.0001, name="b0")# 0.0001, name="b0")

        # def f(x): return tf.cond(tf.less(x, 0.5), lambda: 0, lambda: 1)
        y_pred = tf.add(tf.add(tf.multiply(b1, x1), tf.multiply(b2, x2)), b0)

        loss = tf.reduce_mean(tf.square(y - y_pred))
        # loss = tf.reduce_sum(tf.abs(tf.tanh(y - y_pred)))

        return x1, x2, y, y_pred, loss

def run(X_train, X_test, y_train, y_test):
    x1_batch = np.array(X_train.ratio)
    x2_batch = np.array(X_train.array_length)
    y_batch = np.array(y_train.genotype)

    x1, x2, y, y_pred, loss = linear_regression()

    optimizer = tf.train.GradientDescentOptimizer(0.0000001)#(0.000001)# (0.00000001)
    train_op = optimizer.minimize(loss)

    with tf.Session() as session:
        session.run(tf.global_variables_initializer())
        feed_dict = {
            x1: x1_batch,
            x2: x2_batch,
            y: y_batch
        }

        for i in range(100001):#200001):
            session.run(train_op, feed_dict)
            if i % 1000 == 1:
                print(i, "loss:", loss.eval(feed_dict))

        y_pred = session.run(y_pred, {
            x1: np.array(X_test.ratio),
            x2: np.array(X_test.array_length)
        })

        evaluate(y_test, y_pred)
        # print(b1, b2, b0)

def take_positive(array):
    return np.array([1 if x > 0 else 0 for x in array])


def precision(y, y_pred):
    precision = y.sum() / (y.sum() + y_pred.sum())
    return precision


def specificity(y, y_pred):
    false_negatives = take_positive(np.array((-1 * y_pred) + y, dtype=np.int8)).sum()
    recall = y.sum() / (y.sum() + false_negatives.sum())
    return recall


def recall(y, y_pred):
    false_positives = take_positive(np.array(y_pred - y, dtype=np.int8)).sum()
    specificity = (1 - y).sum() / ((1 - y).sum() + false_positives.sum())
    return specificity


def false_positive_rate(y, y_pred):
    return 1 - recall(y, y_pred)


def false_negative_rate(y, y_pred):
    return 1 - specificity(y, y_pred)


def f_score(y, y_pred):  # best is 1, worst is 0
    return 2 * precision(y, y_pred) * recall(y, y_pred) / (precision(y, y_pred) + recall(y, y_pred))

def fraction_misclassified_00(y, y_pred):
    return take_positive(y_pred - y).sum() / y.sum()

def matthews(y, y_pred):
    TP = recall(y, y_pred)
    TN = specificity(y, y_pred)
    FP = false_positive_rate(y, y_pred)
    FN = false_negative_rate(y, y_pred)
    return (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

def col_normalize(column):  # pd.Series
    mean = column.mean()
    std = column.std()
    new_column = column.apply(lambda x: (x - mean) / std)
    return np.array(new_column)

def round(array):
    return np.array([0 if x < 0.5 else 1 for x in array])

def evaluate(y_test, y_pred):
    '''false_negative = np.array([0 if x >= 0 else 1 for x in (round(y_pred) - y_test)]).sum()
    false_positive = np.array([0 if x >= 0 else 1 for x in (y_test - round(y_pred))]).sum()
    true_negative = np.array([0 if x > 0 else 1 for x in (y_test + round(y_pred))]).sum()
    true_positive = np.array([0 if x <= 1 else 1 for x in (y_test + round(y_pred))]).sum()
    size = len(y_pred)

    # false_negative is what we don't want
    recall = true_positive / (true_positive + false_negative)
    specificity = true_negative / (true_negative + false_positive)
    accuracy = (true_positive + true_negative) / (true_positive + false_positive + true_negative + false_negative)
    precision = true_positive / (true_positive + false_positive)
    matthews_correlation = (true_positive * true_negative - false_positive * false_negative) / np.sqrt(
        (true_positive + false_positive) *
        (true_positive + false_negative) *
        (true_negative + false_positive) *
        (true_negative + false_negative)
    )
    f_measure = 2 * true_positive / (2 * true_positive + false_negative + false_positive)

    print("Total observations: {}\n---\n".format(size))
    print("Recall: {}".format(recall))
    print("Specificity: {}".format(specificity))
    print("Accuracy: {}".format(accuracy))
    print("Precision: {}".format(precision))
    print("Matthews: {}".format(matthews_correlation))
    print("F-measure: {}".format(f_measure))

    print(pd.Series(y_pred).describe())'''

    y_pred = round(y_pred)#np.tanh(y_pred))
    # y_pred = round(y_pred)

    print("---")
    for (yy, zz) in [(precision, "Precision"), (recall, "Recall"),
                     (specificity, "Specificity"), (false_positive_rate, "FPR"),
                     (false_negative_rate, "FNR"), (f_score, "F score"),
                     (fraction_misclassified_00, "Fraction of 0/0's misclassified"), (matthews, "MCC")]:
        print(zz + ": {}".format(yy(y_test, y_pred)))

np.random.seed(100)
df = simulation_set()

df_sub00 = df.loc[df['genotype'] == '0/0']
df_sub01 = df.loc[df['genotype'] == '0/1']
df_sub00 = df_sub00.append(df_sub01)
df_sub00[['genotype']] = df_sub00.genotype.apply(lambda x: 0 if x == '0/0' else 1)
#df_sub00[['genotype']] = df_sub00.genotype.apply(lambda x: 0 if x == '0/0' else 1)

used_columns = ['ratio', 'array_length']
used_labels = ['genotype']

df_sub00[['ratio']] = col_normalize(df_sub00.ratio)
df_sub00[['array_length']] = col_normalize(df_sub00.array_length)

X_train, X_test, y_train, y_test = train_test_split(df_sub00[used_columns], df_sub00[used_labels])

run(X_train, X_test, y_train, np.array(y_test).T.flatten())

