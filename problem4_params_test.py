import pandas as pd
import numpy as np


def simulation_set():
    hg_set = pd.read_csv("/Users/mkorovkin/Desktop/marzd/sim2_148bp.csv")

    hg_set.drop("start", axis=1)
    hg_set.drop("end", axis=1)

    hg_set = hg_set.loc[hg_set["filter"] == "S"]
    # hg_set = hg_set.loc[hg_set["left_flank"] > 0]
    # hg_set = hg_set.loc[hg_set["right_flank"] > 0]
    hg_set = hg_set.loc[hg_set["inside_array"] > 0]

    # Keep flanks static
    lfm = hg_set.left_flank.mean()
    rfm = hg_set.right_flank.mean()
    hg_set["left_flank"] = hg_set.left_flank.apply(lambda x: lfm)
    hg_set["right_flank"] = hg_set.right_flank.apply(lambda x: rfm)

    # Calculate ratios
    hg_set["ratio"] = hg_set.inside_array / (hg_set.left_flank + hg_set.right_flank)

    # prepare statistical data
    sm = pd.read_csv("/Users/mkorovkin/Desktop/output_statistics_00_new.csv")
    xs = np.arange(148 * 3, 148 * 40, 148)
    ms = sm.mA148

    # analyze total_set
    hg_set['inside_length_ratio'] = hg_set.inside_array / hg_set.array_length

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

    ven = {
        "43466": 1,
        "1": 1,
        "0/1": 1,
        "0/-1": 1,
        "0/0": -1,
    }

    hg_set['label'] = hg_set.genotype.apply(lambda x: ven[x])

    return hg_set

def precision(y, y_pred):
    true_positives = np.array([1 if (ya == 1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)])
    false_positives = np.array([1 if (ya == -1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)])
    return true_positives.sum() / (true_positives.sum() + false_positives.sum())

def number_fp(y, y_pred):
    return np.array([1 if (ya == -1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)]).sum()

def number_tp(y, y_pred):
    return np.array([1 if (ya == 1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)]).sum()

def number_fn(y, y_pred):
    return np.array([1 if (ya == 1 and yp == -1) else 0 for (ya, yp) in zip(y, y_pred)]).sum()

def number_tn(y, y_pred):
    return np.array([1 if (ya == -1 and yp == -1) else 0 for (ya, yp) in zip(y, y_pred)]).sum()

def recall(y, y_pred):
    # false_negatives = take_positive(np.array((-1 * y_pred) + y, dtype=np.int8)).sum()
    # recall = y.sum() / (y.sum() + false_negatives.sum())
    true_positives = np.array([1 if (ya == 1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)])
    false_negatives = np.array([1 if (ya == 1 and yp == -1) else 0 for (ya, yp) in zip(y, y_pred)])
    return true_positives.sum() / (true_positives.sum() + false_negatives.sum())


def specificity(y, y_pred):
    # false_positives = take_positive(np.array(y_pred - y, dtype=np.int8)).sum()
    # specificity = (1 - y).sum() / ((1 - y).sum() + false_positives.sum())
    # return specificity
    true_negatives = np.array([1 if (ya == -1 and yp == -1) else 0 for (ya, yp) in zip(y, y_pred)])
    false_positives = np.array([1 if (ya == -1 and yp == 1) else 0 for (ya, yp) in zip(y, y_pred)])
    return true_negatives.sum() / (true_negatives.sum() + false_positives.sum())

def false_positive_rate(y, y_pred):
    return 1 - specificity(y, y_pred)

def false_negative_rate(y, y_pred):
    return 1 - recall(y, y_pred)

def f_score(y, y_pred):  # best is 1, worst is 0
    return 2 * precision(y, y_pred) * recall(y, y_pred) / (precision(y, y_pred) + recall(y, y_pred))

def matthews(y, y_pred):
    TP = recall(y, y_pred)
    TN = specificity(y, y_pred)
    FP = false_positive_rate(y, y_pred)
    FN = false_negative_rate(y, y_pred)
    return (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

def output_results(y_actual, preds):
    print("---{}---".format("Results"))
    for (yy, zz) in [(precision, "Precision"), (recall, "Recall"),
                     (specificity, "Specificity"), (false_positive_rate, "FPR"),
                     (false_negative_rate, "FNR"), (f_score, "F score"),
                     (number_tp, "True positive count",), (number_fp, "False positive count"),
                     (number_tn, "True negative count"), (matthews, "MCC")]:
        print(zz + ": {}".format(yy(y_actual, preds)))


def predict_from_tf_SVCNeg(x_y_zip, params):
    slope = (-params[0] / params[1])
    intercept = (params[2] / params[1])

    class_result = [1 if (x * slope + intercept > y) else -1 for (x, y) in x_y_zip]
    return np.array(class_result)


def predict_from_tf_SVCPos(x_y_zip, params):
    slope = (-params[0] / params[1])
    intercept = (params[2] / params[1])

    class_result = [1 if (x * slope + intercept < y) else -1 for (x, y) in x_y_zip]
    return np.array(class_result)

dataset = simulation_set()

x_columns = ['array_length', 'ratio']
y_columns = ['label']

dataset_X = dataset[x_columns]
dataset_y = dataset[y_columns]

params_neg = [-0.0254627823167, 4.308455495333334, -3.793559205499995]
params_pos = [-0.0338727795133, 4.294555499966664, -4.498059203999996]

preds_neg = predict_from_tf_SVCNeg(zip(np.array(dataset_X.array_length), np.array(dataset_X.ratio)), params_neg)
preds_pos = predict_from_tf_SVCPos(zip(np.array(dataset_X.array_length), np.array(dataset_X.ratio)), params_pos)
classification_list = [1 if (x == 1 or y == 1) else -1 for (x, y) in zip(preds_neg, preds_pos)]

output_results(np.array(dataset_y.label), classification_list)