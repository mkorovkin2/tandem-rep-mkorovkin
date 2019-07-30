import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from scipy.stats import linregress
import tensorflow as tf

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

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

#################################################### STATISTICS ########################################################

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

###################################################### SVM IDEA ########################################################

def run_linear_SVM(X_train, X_test, y_train, y_test, X_total, y_total, batch_size=96, alpha=0.01, C=1.0, d=0.01, d_penalty=1.0, bound_mod=1.0, lr=0.001, fp_penalty=1., epochs=100000):
    # np.max(len(X_train) // 100, 96)

    x_data = tf.placeholder(shape=(None, 2), dtype=tf.float32)
    y_label = tf.placeholder(shape=(None, 1), dtype=tf.float32)

    # [[a0], [a1]] = [[2.8582766], [-0.020608412]]
    # [[b]] = [[-2.540259]]
    #[[2.8826115], [-0.02194276]] /// [[2.8582766], [-0.022473007]]
    #b: [[-2.5355039]] /// [[-2.540259]]

    #>> A: [[ 2.8860528 ]
 #[-0.02145437]]
#>> b: [[-2.5345695]]

    #[[4.2944555 + a0mod], [-0.03241278 + a1mod]]
    #[[b]] = [[-3.9100592 + bmod]]

    A = tf.Variable([[4.2944555], [-0.03241278]])#[[4.2805715], [-0.032344025]])# [[2.8860528 ], [-0.02145437]])
    #[[2.8826115], [-0.02194276]])#[[2.8582766], [-0.020608412]])#tf.divide(tf.random_normal(shape=[2, 1]), 10))
    b = tf.Variable([[-3.9100592 - 0.6]])#[[-3.8925476]])# [[-2.5345695]])
    #[[-2.5355039]])#[[-2.540259]])#tf.divide(tf.random_normal(shape=[1, 1]), 10))



    model_output = tf.subtract(
        tf.matmul(
            x_data, A
        ), b
    )
    l2_norm = tf.reduce_sum(tf.square(A))

    alpha = tf.constant([alpha])
    C = tf.constant([C])
    d = tf.constant([d])
    bound_mod = tf.constant([bound_mod])

    # classification_term = tf.reduce_sum(
    #    tf.maximum(
    #         0., tf.subtract(
    #             1., tf.multiply(
    #                 model_output, y_label
    #             )
    #         )
    #     )
    # )

    error_term = tf.reduce_sum(
        tf.maximum(
            0., tf.subtract(
                1.0, tf.multiply(
                    y_label, model_output
                )
            )
        )
    )

    error_term_penalty = tf.reduce_sum(
        tf.maximum(
            0., tf.subtract(
                1.0, tf.multiply(
                    tf.subtract(y_label, d), model_output
                )
            )
        )
    )

    '''
    error_term_positive_deter = tf.reduce_sum(
        tf.maximum(
            tf.maximum(d, 0), tf.subtract(
                tf.multiply(1., bound_mod), tf.multiply(
                    model_output, tf.subtract(y_label, d)
                )
            )
        )
    )
    
    error_term_positive_deter = tf.reduce_sum(
        tf.maximum(
            tf.maximum(d, 0), tf.subtract(
                tf.multiply(1., bound_mod), tf.multiply(
                    model_output, tf.subtract(y_label, d)
                )
            )
        )
    )

    error_residuals = tf.reduce_sum(
        tf.subtract(
            1.0, tf.add(
                error_term, tf.multiply(
                    model_output, y_label
                )
            )
        )
    )
    '''

    prediction = tf.sign(model_output)
    accuracy = tf.reduce_mean(
        tf.cast(
            tf.equal(
                prediction, y_label
            ), tf.float32
        )
    )

    true_positives = tf.reduce_sum(tf.cast(tf.equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.float32))
    false_positives = tf.reduce_sum(tf.cast(tf.equal(tf.not_equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.equal(prediction, 1.0)), tf.float32))
    true_negatives = tf.reduce_sum(tf.cast(tf.equal(tf.equal(y_label, -1.0), tf.equal(prediction, -1.0)), tf.float32))
    false_negatives = tf.reduce_sum(tf.cast(tf.equal(tf.not_equal(tf.equal(y_label, -1.0), tf.equal(prediction, -1.0)), tf.equal(prediction, -1.0)),tf.float32))

    precision_result = tf.divide(true_positives, tf.add(true_positives, false_positives))
    recall_result = tf.divide(true_positives, tf.add(true_positives, false_negatives))
    specificity_result = tf.divide(true_negatives, tf.add(false_positives, true_negatives))

    '''
    precision_result = tf.divide(tf.reduce_sum(tf.cast(tf.equal(tf.not_equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.equal(prediction, 1)), tf.float32)),
        tf.add(
            tf.reduce_sum(tf.cast(tf.equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.float32)),
            tf.reduce_sum(tf.cast(tf.equal(tf.not_equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.equal(prediction, 1)), tf.float32))
        )
    )

    recall_result = tf.divide(tf.reduce_sum(
        tf.cast(tf.equal(tf.equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.equal(prediction, 1)),
                tf.float32)),
                          tf.add(
                              tf.reduce_sum(
                                  tf.cast(tf.equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)), tf.float32)),
                              tf.reduce_sum(tf.cast(
                                  tf.equal(tf.not_equal(tf.equal(y_label, 1.0), tf.equal(prediction, 1.0)),
                                           tf.equal(prediction, 1)), tf.float32))
                          )
                          )

    fp_result = tf.subtract(1., recall_result)

    fn_result = tf.divide(tf.reduce_sum(
        tf.cast(tf.equal(tf.not_equal(tf.equal(y_label, -1.0), tf.equal(prediction, -1.0)), tf.equal(prediction, -1.0)),
                tf.float32)),
                          tf.add(
                              tf.reduce_sum(
                                  tf.cast(tf.equal(tf.equal(y_label, -1.0), tf.equal(prediction, -1.0)), tf.float32)),
                              tf.reduce_sum(tf.cast(
                                  tf.equal(tf.not_equal(tf.equal(y_label, -1.0), tf.equal(prediction, -1.0)),
                                           tf.equal(prediction, -1.0)), tf.float32))
                          )
                          )
    '''

    loss = tf.add(tf.maximum(0.0, tf.multiply(C, error_term_penalty)), l2_norm)

    '''
    loss = tf.add(
        #tf.add(
        tf.multiply(
            C, error_term
        ), tf.add(
            l2_norm, #tf.add(l2_norm, tf.multiply(
            #    error_term_positive_deter, d_penalty
            #)),
            tf.multiply(
                alpha, error_residuals
            )
        )
    )#, tf.multiply(fp_result, fp_penalty))
    '''

    optimizer = tf.train.AdamOptimizer(lr)#tf.train.GradientDescentOptimizer(lr)# AdamOptimizer(lr)# 0.01)# 0.001) # tf.train.GradientDescentOptimizer(0.00000001)
    training_step = optimizer.minimize(loss)

    init = tf.global_variables_initializer()
    sess.run(init)

    loss_vec = []
    train_accuracy = []
    test_accuracy = []

    for i in range(epochs):
        #if batch_size <= 0:
        #    rand_x = X_train  # [rand_index]
        #    rand_y = y_train  # [rand_index]
        #else:
        rand_index = np.random.choice(len(X_train), size=batch_size)
        rand_x = X_train[rand_index]
        rand_y = y_train[rand_index]


        sess.run(training_step, feed_dict={
            x_data: rand_x,
            y_label: rand_y
        })

        current_loss = sess.run(loss, feed_dict={
            x_data: rand_x,
            y_label: rand_y
        })
        loss_vec.append(current_loss)

        precision_train = sess.run(precision_result, feed_dict={
            x_data: rand_x,
            y_label: rand_y
        })

        recall_train = sess.run(recall_result, feed_dict={
            x_data: rand_x,
            y_label: rand_y
        })

        specificity_train = sess.run(specificity_result, feed_dict={
            x_data: rand_x,
            y_label: rand_y
        })

        # print("Loss: {}".format(current_loss))

        current_train_accuracy = sess.run(accuracy, feed_dict={
            x_data: X_train,
            y_label: y_train
        })
        train_accuracy.append(current_train_accuracy)

        '''
        current_test_accuracy = sess.run(accuracy, feed_dict={
            x_data: X_test,
            y_label: y_test
        })
        test_accuracy.append(current_test_accuracy)

        false_positive_test = sess.run(fp_result, feed_dict={
            x_data: X_test,
            y_label: y_test
        })
        '''

        if (i + 1) % 1000 == 0:
            '''
            print('Step #{}\n  Loss = {}\n  Train accuracy = {}\n  False positive rate = {}'.format(
                str(i + 1),
                current_loss,
                current_train_accuracy,
                false_positive_train
            ))
            '''
            print("Step #{} | Recall = {} | Precision = {} | FPR {} | loss = {} | accuracy = {}".format(
                str(i + 1),
                recall_train,
                precision_train,
                1 - specificity_train,
                current_loss,
                current_train_accuracy))

            if recall_train > 0.78:
                print(">>> A: {}".format(np.array(sess.run(A)).squeeze()))
                print(">>> b: {}".format(sess.run(b)))

    [[a0], [a1]] = sess.run(A)
    [[b]] = sess.run(b)
    print("A:", [[a0], [a1]])
    print("b:", [[b]])
    print("Test recall:", sess.run(recall_result, feed_dict={
        X_test,
        y_test
    }))
    print("Test FNR:", sess.run(false_negatives, feed_dict={
        X_test,
        y_test
    }))
    print("Test acc:", sess.run(accuracy, feed_dict={
        X_test,
        y_test
    }))

    svm_slope = -a1 / a0
    y_intercept = b / a0
    x_total_values = [d[1] for d in X_total]
    best_fit = []
    for i in x_total_values:
        best_fit.append(svm_slope * i + y_intercept)

    X_set1 = [d[1] for i, d in enumerate(X_total) if y_total[i] == 1]
    y_set1 = [d[0] for i, d in enumerate(X_total) if y_total[i] == 1]
    X_set0 = [d[1] for i, d in enumerate(X_total) if y_total[i] == -1]
    y_set0 = [d[0] for i, d in enumerate(X_total) if y_total[i] == -1]

    # Plot everything
    plt.plot(X_set1, y_set1, 'o', label='0/1 Set')
    plt.plot(X_set0, y_set0, 'x', label='0/0 Set')
    plt.plot(x_total_values, best_fit, 'r-', label='Linear Separator', linewidth=1)
    plt.legend(loc='lower right')
    plt.title('Array Length vs. Ratio')
    plt.xlabel('Array Length')
    plt.ylabel('Ratio')
    plt.xlim([148, 2500])
    plt.ylim([0, 19])
    plt.show()

    # Plot train/test accuracies
    plt.plot(train_accuracy, 'k-', label='Training Accuracy')
    plt.plot(test_accuracy, 'r--', label='Test Accuracy')
    plt.title('Train and Test Set Accuracies')
    plt.xlabel('Generation')
    plt.ylabel('Accuracy')
    plt.legend(loc='lower right')
    plt.show()

    # Plot loss over time
    plt.plot(loss_vec, 'k-')
    plt.title('Loss per Generation')
    plt.xlabel('Generation')
    plt.ylabel('Loss')
    plt.show()

#################################################### SCRIPT CODE #######################################################

np.random.seed(100)
tf.set_random_seed(41)
sess = tf.Session()

df = simulation_set()

df_sub00 = df.loc[df['genotype'] == '0/0']
df_sub01 = df.loc[df['genotype'] == '0/1']
df_sub00 = df_sub00.append(df_sub01)
df_sub00[['genotype']] = df_sub00.genotype.apply(lambda x: -1 if x == '0/0' else 1)

used_columns = ['ratio', 'array_length']
used_labels = ['genotype']

# df_sub00[['ratio']] = col_normalize(df_sub00.ratio)
# df_sub00[['array_length']] = col_normalize(df_sub00.array_length)

X_train, X_test, y_train, y_test = train_test_split(df_sub00[used_columns], df_sub00[used_labels])

run_linear_SVM(X_train.values,
               X_test.values,
               y_train.values,
               y_test.values,
               df_sub00[used_columns].values,
               df_sub00[used_labels].values,
               batch_size=54,
               # alpha=0.01,
               C=10.0,
               d=0.3, #-0.7, # -0.7,
               # d_penalty=0.001,#1.0, # 100.0,
               # bound_mod=1.0,
               lr=0.001,#large batch0.01,#ProximalAdagrad - 0.001,# THIS IS FOR THE BATCH 256 - 0.000000001, # 0.0000001, # 0.0000000001,
               # fp_penalty=100.0,
               epochs=100000)

# problems coming from