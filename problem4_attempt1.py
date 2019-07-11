import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from scipy.stats import linregress
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import BaggingRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor

pd.set_option('display.max_columns', 500)

def use_ref_set():
    ref_set = pd.read_csv("/Users/mkorovkin/Desktop/marzd/refset_full.csv")
    hg_set = pd.read_csv("/Users/mkorovkin/Desktop/marzd/sim2_148bp.csv")
    # pd.read_csv("/Users/mkorovkin/Desktop/marzd/HG001counts.csv") # read length = 148

    # Merge the ref_set and the hg_set on TRID
    ref_set = ref_set.loc[ref_set["indistinguishable"] == "S"]
    ref_set[["TRID"]] = ref_set.TRID.apply(lambda x: str(x) + "_S")

    hg_set = hg_set.loc[hg_set["TRID"].str.contains("_S")]
    hg_set = hg_set.loc[hg_set["left_flank"] > 0]
    hg_set = hg_set.loc[hg_set["right_flank"] > 0]
    hg_set = hg_set.loc[hg_set["inside_array"] > 0]

    hg_set["ratio"] = hg_set.inside_array / (hg_set.left_flank + hg_set.right_flank)

    ref_set_columns = ["TRID", "pattern_size", "copy_number"]
    ref_set = ref_set[ref_set_columns]

    # total_set is now the merged set
    total_set = pd.merge(hg_set, ref_set, on="TRID")

    # prepare statistical data
    sm = pd.read_csv("/Users/mkorovkin/Desktop/output_statistics_00_new.csv")
    xs = np.arange(148 * 3, 148 * 40, 148)
    ms = sm.mA148

    slope1, intercept1, _, _, _ = linregress(xs, ms)

    # analyze total_set
    total_set['inside_length_ratio'] = total_set.inside_array / total_set.array_length
    total_set['error'] = (total_set.array_length.apply(lambda x: slope1 * x + intercept1) - total_set.ratio) / total_set.ratio

    outliers = total_set.loc[total_set['inside_length_ratio'] < 0.1]

    print(
        total_set.loc[total_set['inside_length_ratio'] < 0.1]
    )
    print("---")
    print(
        total_set.loc[total_set['inside_length_ratio'] > 6]
    )

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
    # print("flanks", lfm, rfm)
    hg_set["left_flank"] = hg_set.left_flank.apply(lambda x: lfm)
    hg_set["right_flank"] = hg_set.right_flank.apply(lambda x: rfm)

    # Calculate ratios
    hg_set["ratio"] = hg_set.inside_array / (hg_set.left_flank + hg_set.right_flank)

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

    '''
    return
    outliers = hg_set.loc[hg_set['inside_length_ratio'] < 0.1]

    print(
        hg_set.loc[hg_set['inside_length_ratio'] < 0.1]
    )
    print("---")
    plt.scatter(outliers.inside_length_ratio, outliers.error)
    plt.show()
    '''

def score_direct(preds, correct):
    cor = 0
    plist = []
    ylist = []
    for (p, y) in zip(preds, correct):
        if (p == y):
            cor += 1
        else:
            plist.append(p)
            ylist.append(y)
    # print("plist\n", pd.Series(plist).value_counts())
    # print("ylist\n", pd.Series(ylist).value_counts())
    return cor / len(correct) if len(correct) > 0 else 0

def score_mae(preds, correct):
    err = (preds - correct).mean()
    return np.array(abs(preds - correct - err)).mean()

def score_mae_sim(preds, correct):
    preds2 = np.array(np.round(preds, decimals=0), dtype=np.int8)
    err = (preds2 - correct).mean()
    preds = np.array(np.round(preds2 - err, decimals=0), dtype=np.int8)
    return score_direct(preds, correct)

def score_ml(preds, correct):
    preds = close_round(preds, correct)
    return score_direct(preds, correct)

def score_individually(preds, correct):
    preds = np.array(np.round(preds, decimals=0), dtype=np.int8)

    print("-2 ::", np.round(score_direct(preds[correct == -2], correct[correct == -2]), decimals=4))
    print("-1 ::", np.round(score_direct(preds[correct == -1], correct[correct == -1]), decimals=4))
    print("0  ::", np.round(score_direct(preds[correct == 0], correct[correct == 0]), decimals=4))
    print("1  ::", np.round(score_direct(preds[correct == 1], correct[correct == 1]), decimals=4))
    print("2  ::", np.round(score_direct(preds[correct == 2], correct[correct == 2]), decimals=4))

    # print(preds[correct == 2], correct[correct == 2])

# ## # ## # # # # ## # ## # # ## ## # ## # # ## ## # ## # # ## ## # ## # # ## ## # ## # # ## ## # ## # # ## ## # ## # #

def close_round(preds, correct, epochs=100, learning_rate=0.001, decimals=0, modify=True, score_func=score_mae):
    bounds_dep = {
        "-2top": -1.4,
        # "-1bot": -1.4,
        # "-1top": -0.4,
        "0bot": -0.4,
        # "0top": 0.4,
        "1bot": 0.4,
        # "1top": 1.4,
        "2bot": 1.4,
    }
    bounds = {
        -2: -1.4,
        -1: -0.4,
        0: 0.4,
        1: 1.4
    }

    if modify:
        for key in bounds.keys():
            bounds[key] += np.random.uniform(-0.1, 0.1)

    for i in range(epochs):
        score, preds_results, misclass = close_round_helper(preds, bounds, correct, score_func)
        bounds, flag = fix_misclassify(score, preds_results, misclass, bounds, learning_rate, score_func, correct)

        print("Epoch " + str(i) + ":", score)

        if flag:
            print("\nComplete with early stopping.\n", score)
            break

    return bounds

def close_round_helper(preds, bounds, correct, score_func):
    preds2 = preds.copy()

    mask0 = np.logical_and(preds < bounds[0], preds > bounds[-1])
    mask1 = np.logical_and(preds <= bounds[1], preds >= bounds[0])
    maskn1 = np.logical_and(preds <= bounds[-1], preds >= bounds[-2])
    preds2[mask0] = 0
    preds2[mask1] = 1
    preds2[preds > bounds[1]] = 2
    preds2[maskn1] = -1
    preds2[preds < bounds[-2]] = -2

    return score_func(preds2, correct), preds2, misclassified(preds, correct)

def fix_misclassify(score, pred_results, misclass, bounds, learning_rate, score_func, correct):
    # find top misclassified
    sorted_misclass = sorted((value, key) for (key, value) in misclass.items())

    top_list = [
        sorted_misclass[0][1], sorted_misclass[1][1], sorted_misclass[2][1]
    ]
    mod_list = [-1, 0, 1]

    running_min = score
    id = (1, 1, 1)
    adjust_bounds = adjacent(sorted(top_list))
    # print("Adjust", adjust_bounds)

    new_mod_dict = {}
    for i in mod_list:
        for j in mod_list:
            for k in mod_list:
                bounds[adjust_bounds[0]] += mod_list[i] * learning_rate
                if len(adjust_bounds) > 1:
                    bounds[adjust_bounds[1]] += mod_list[j] * learning_rate
                if len(adjust_bounds) > 2:
                    bounds[adjust_bounds[2]] += mod_list[k] * learning_rate

                sc1, ps2, msf = close_round_helper(pred_results, bounds, correct, score_func)
                if sc1 < running_min:
                    running_min = sc1
                    id = (i, j, k)

                new_mod_dict[i] = {j:
                                       {k: sc1}
                                   }
                bounds[adjust_bounds[0]] -= mod_list[i] * learning_rate
                if len(adjust_bounds) > 1:
                    bounds[adjust_bounds[1]] -= mod_list[j] * learning_rate
                if len(adjust_bounds) > 2:
                    bounds[adjust_bounds[2]] -= mod_list[k] * learning_rate

    bounds[adjust_bounds[0]] += mod_list[id[0]] * learning_rate
    if len(adjust_bounds) > 1:
        bounds[adjust_bounds[1]] += mod_list[id[1]] * learning_rate
    if len(adjust_bounds) > 2:
        bounds[adjust_bounds[2]] += mod_list[id[2]] * learning_rate

    flag = False
    if id[0] == 1 and id[1] == 1 and id[2] == 1:
        flag = True

    return bounds, flag

def adjacent(tlist):
    ret_list = []
    last_added = False
    '''for i in range(len(tlist) - 1):
        if abs(tlist[i] - tlist[i + 1]) <= 1.1:
            if not last_added:
                ret_list.append(tlist[i])
            else:
                ret_list.append(tlist[i + 1])
            last_added = False
        else:
            if not last_added:
                ret_list.append(tlist[i])
                ret_list.append(tlist[i + 1])
                last_added = True
            else:
                ret_list.append(tlist[i + 1])'''

    for i in tlist:
        if i < 2:
            ret_list.append(i)
    return ret_list

def misclassified(preds, correct):
    return {
        0: score_direct(correct[correct == 0], preds[correct == 0]),
        1: score_direct(correct[correct == 1], preds[correct == 1]),
        2: score_direct(correct[correct == 2], preds[correct == 2]),
        -1: score_direct(correct[correct == -1], preds[correct == -1]),
        -2: score_direct(correct[correct == -2], preds[correct == -2])
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dataset = simulation_set()
print(list(dataset))

x_columns = ['array_length', 'pattern_size', 'array_length', 'copy_number', 'error', 'ratio']
y_columns = ['label']

dataset_X = dataset[x_columns]
dataset_y = dataset[y_columns]

train_X, test_X, train_y, test_y = train_test_split(dataset_X, dataset_y, test_size=0.2)

print(train_X.shape)

rfr = RandomForestRegressor(random_state=100)
svc = SVC(random_state=100)
gpr = GaussianProcessRegressor(random_state=100)
knn = KNeighborsClassifier()
bcl = BaggingClassifier(random_state=100)
gbc = GradientBoostingClassifier(random_state=100)
gbr = GradientBoostingRegressor(random_state=100)

for (reg, lab) in (rfr, "RandomForestRegressor"),\
                  (svc, "SVC"),\
                  (gpr, "GaussianProcessRegressor"),\
                  (knn, "KNeighborsClassifier"),\
                  (bcl, "BaggingClassifier"),\
                  (gbc, "GradientBoostingClassifier"),\
                  (gbr, "GradientBoostingRegressor"):
    classifier = reg.fit(train_X, train_y.label)
    preds = classifier.predict(test_X)

    print("---{}---".format(lab))
    print("{} direct accuracy score: {}".format(lab,
                                                np.round(score_direct(preds, test_y.label), decimals=4)
                                                ))
    print("{} MAE score: {}".format(lab,
                                    np.round(score_mae(preds, test_y.label), decimals=4)
                                    ))
    print("{} adjusted MAE score: {}".format(lab,
                                             np.round(score_mae_sim(preds, np.array(test_y.label)), decimals=4)
                                             ))
    print("Breakdown of individual accuracies for {}".format(lab))
    score_individually(preds, np.array(test_y.label))
    print("\n")

    # print("{} bounds: {}".format(lab, close_round(preds, np.array(test_y.label), learning_rate=0.0001)))
    # print(adjacent([-1, 1, 2]))
    # print(adjacent([0, 1, 2]))