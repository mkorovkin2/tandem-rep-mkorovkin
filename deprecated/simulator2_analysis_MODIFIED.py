import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import linregress

def sim_reg_aux(A, R, F, t, sim_times, C):
    T = A + 2 * F
    N = np.int32(np.ceil(C * T / R))

    hit_count_list = list()
    out_count_list = list()

    for _ in range(sim_times):
        hit_count = 0
        out_count = 0
        starting_points = list()

        A_lower = F
        A_upper = (A_lower - 0.25 * F) + A - 0.75 * F - R
        F_lower = 0
        F_lower_bound = F * 0.5
        F_upper = A_upper
        F_upper_bound = F_upper + F * 0.5

        for _ in range(np.int32(N / 2)):
            xs = np.floor(rand.uniform(0, A_upper + np.int32(0.5 * F)))
            starting_points.append(xs)
            starting_points.append(xs + t)

        for i in range(len(starting_points)):
            x = starting_points[i]
            if (F_lower <= x and F_lower_bound > x) or (F_upper < x and F_upper_bound >= x):
                hit_count += 1
            elif (A_lower <= x and A_upper >= x):
                out_count += 1

        hit_count_list.append(hit_count)
        out_count_list.append(out_count)

    sum_hit = 0
    sum_out = 0
    ratio_list = list()

    for i in range(sim_times):
        sum_hit += hit_count_list[i]
        sum_out += out_count_list[i]
        if hit_count_list[i] > 0:
            ratio_list.append(out_count_list[i] / hit_count_list[i])

    std_return = np.round(np.std(ratio_list), decimals=2)
    mean_return = np.round(np.mean(ratio_list), decimals=2)
    print("A=" + str(A) + ", C=" + str(C) + ", mean=" + str(mean_return) + ", std=" + str(std_return))
    return mean_return, std_return

def find_intersections(mAlist, sAlist, Alist):
    adjacency_matrix = np.zeros(shape=(len(mAlist), len(mAlist)))
    for mi in range(len(mAlist)):
        std = sAlist[mi]
        mean = mAlist[mi]
        for smi in range(len(mAlist)):
            substd = sAlist[smi]
            submean = mAlist[smi]
            if (mean + std < submean - substd) or (mean - std > submean + substd):
                adjacency_matrix[mi][smi] = Alist[smi]
            else:
                adjacency_matrix[mi][smi] = 0
    return adjacency_matrix

def find_variance(kwargs, A_v, R_v, F_v):
    to_find = "variance"
    followup_variance = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip()
    if "in" in followup_variance:
        to_find = "in"
        followup_in = followup_variance[(followup_variance.index(to_find) + len(to_find) + 1):].strip()
        if "A" in followup_in:
            to_find = "A"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            A_v = (True, np.float32(followup[0]))
        elif "R" in followup_in:
            to_find = "R"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            R_v = (True, np.float32(followup[0]))
        elif "F" in followup_in:
            to_find = "F"
            followup = followup_in[(followup_in.index(to_find) + len(to_find) + 1):].strip().split(" ")
            F_v = (True, np.float32(followup[0]))
    if "variance" in followup_variance:
        return find_variance(followup_variance, A_v, R_v, F_v)
    else:
        return A_v, R_v, F_v

def sim_adv(A, R, F, sim_times, C, fragment_imp):
    T = A + 2 * F
    N = np.int32(np.ceil(C * T / R))

    hit_count_list = list()
    out_count_list = list()

    for _ in range(sim_times):
        hit_count = 0
        out_count = 0
        starting_points = list()

        random_normal_values = np.random.normal(fragment_imp[0], fragment_imp[1], np.int32(N / 2))
        for index in range(len(random_normal_values)):
            if np.floor(random_normal_values[index]) > 0:
                random_index1 = np.ceil(np.floor(rand.uniform(0, A + np.int32(0.5 * F))))
                random_index2 = np.floor(random_index1 + random_normal_values[index] - R)
                starting_points.append(random_index1)
                starting_points.append(random_index2)

        A_lower = F
        A_upper = (A_lower - 0.25 * F) + A - 0.75 * F - R
        F_lower = 0
        F_lower_bound = F * 0.5
        F_upper = A_upper
        F_upper_bound = F_upper + F * 0.5

        for i in range(len(random_normal_values)):
            x = starting_points[i]
            if (F_lower <= x and F_lower_bound > x) or (F_upper < x and F_upper_bound >= x):
                hit_count += 1
            elif (A_lower <= x and A_upper >= x):
                out_count += 1

        hit_count_list.append(hit_count)
        out_count_list.append(out_count)

    sum_hit = 0
    sum_out = 0
    ratio_list = list()

    for i in range(sim_times):
        sum_hit += hit_count_list[i]
        sum_out += out_count_list[i]
        if hit_count_list[i] > 0:
            ratio_list.append(out_count_list[i] / hit_count_list[i])

    std_return = np.round(np.std(ratio_list), decimals=2)
    mean_return = np.round(np.mean(ratio_list), decimals=2)
    print("A=" + str(A) + ", C=" + str(C) + ", mean=" + str(mean_return) + ", std=" + str(std_return))
    return mean_return, std_return

def runsim(Clist, READ_LENGTH, fragment_variance, sim_times, test_list):
    mlist = list()
    slist = list()
    tlist = list()

    mdict = {}
    Alist = np.arange(READ_LENGTH * 2, READ_LENGTH * 8, READ_LENGTH)

    for C in Clist:
        mdict[C] = {}
        for A in Alist:
            mdict[C][A] = (0, 0)

    ordered_mean = list()
    ordered_ci = list()
    ordered_t = list()
    iter = 0

    for t in test_list:
        for A in Alist:
            cmean = 0
            cstd = 0
            for C in Clist:
                mean, std = sim_adv(A, READ_LENGTH, READ_LENGTH, sim_times, C, (t, 150))
                mdict[C][A] = (mean, std)
                cmean += mean
                cstd += std

            cmean = cmean / len(Clist)
            cstd = cstd / len(Clist)
            ordered_mean.append(cmean)
            ordered_ci.append(cstd)

            ordered_t.append(A)
            print("> Progress: " + str(100 * np.round(iter / (len(test_list) * len(Alist)), decimals=2)) + "%")
            iter += 1

    color_wheel = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    for index in range(len(Alist)):
        plt.errorbar(test_list, ordered_mean[index::len(Alist)], capsize=4, yerr=ordered_ci[index::len(Alist)], color=color_wheel[index])

    plt.title("READ_LENGTH=" + str(READ_LENGTH) + " | " + str(sim_times) + " iterations per simulation" + (" | COVERAGE=" + str(Clist[0]) if len(Clist) < 2 else ""))
    plt.ylabel("X/Y ratio")
    plt.xlabel("fragment length | std.dev=" + str(fragment_variance))
    patches1 = [patches.Patch(color=color_wheel[Ai], label=("Array length: " + str(Alist[Ai]))) for Ai in range(len(Alist))]
    plt.legend(handles=patches1)
    plt.xlim(np.int32(min(test_list) * 24 / 25), np.int32(max(test_list) * 5 / 4))
    plt.show()

def sim_regression_with_randomized_fragment(Clist, fragment_variance, fragment_length_mean, sim_times, sim_times_rand, array_length, read_length, flank_size, plot):
    #test_list = np.array(np.random.normal(fragment_length_mean, fragment_variance, sim_times_rand), dtype=np.int32)
    test_list = np.array(np.random.uniform(read_length, fragment_length_mean * 2 + (fragment_length_mean - read_length), sim_times_rand), dtype=np.int32)
    ordered_mean = list()
    ordered_ci = list()
    ordered_t = list()

    iter = 0
    for t in test_list:
        cmean = 0
        cstd = 0
        for C in Clist:
            mean, std = sim_reg_aux(array_length, read_length, flank_size, t, sim_times, C)
            cmean += mean
            cstd += std

        cmean = cmean / len(Clist)
        cstd = cstd / len(Clist)
        ordered_mean.append(cmean)
        ordered_ci.append(cstd * 1.96 / np.sqrt(sim_times))
        ordered_t.append(t)

        print("> Progress: " + str(100 * np.round(iter / (len(test_list)), decimals=2)) + "%")
        iter += 1

    #color_wheel = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    slope1, intercept1, r_value1, p_value1, std_err1 = linregress(ordered_t, ordered_mean)
    print("Beta1:", slope1)
    print("Beta2:", intercept1)
    print("R-squared:", r_value1 ** 2)
    print("P-value:", p_value1)
    print("Std.err:", std_err1)
    if plot:
        plt.errorbar(ordered_t, ordered_mean, yerr=ordered_ci, elinewidth=1, capsize=1, color='b', fmt='o')
        maxx = np.int32(max(ordered_t) * 0.9)
        minx = np.int32(min(ordered_t) * 1.1)
        plt.title("array_length=" + str(array_length) + " | read_length=" + str(read_length) + " | iterations=" + str(sim_times * sim_times_rand))
        plt.xlabel("mean fragment length = " + str(fragment_length_mean) + " | std.dev = " + str(fragment_variance))
        plt.ylabel("X/Y ratio of sample mean")
        x_vals_reg = np.arange(minx, maxx, 1)
        plt.plot(x_vals_reg, [(slope1 * x + intercept1) for x in x_vals_reg], color='r')
        plt.show()

    cimean = np.array(ordered_ci).mean() * 2
    (intercept1 - cimean)

Clist = [70, 355]
READ_LENGTH = 148#100#250
fragment_variance = 150
sim_times = 20
fragment_length_list = [400, 550, 700]

fragment_length_mean = 550
array_length_sim_reg = 1500
beta = sim_regression_with_randomized_fragment(
    Clist, fragment_variance,
    fragment_length_mean, sim_times,
    100,#np.int32(sim_times / 4) if sim_times > 60 else 30,
    array_length_sim_reg, READ_LENGTH,
    np.int32(READ_LENGTH * 0.75), True
)

exit(0)

for C in Clist:
    runsim([C],
           READ_LENGTH,
           fragment_variance,
           sim_times,
           fragment_length_list
           )

exit(0)

# Fragment size from left end to right end; read lengths [100, 150, 250]; separation of coverages; increase number of simulations to 100;
#       allow fragment length to be smaller than 2R; fragment lengths 400-700 (550~150);
# Graph of paired reads with error bars; -> add regression line; regression on 95% confidence intervals
# Describe simulations in text on a document
# loss: losing data from re-mapping reads (chance of loss?)

# ultimately: how different do ratios have to be so that you can tell that array length changed; can you tell the length of the array (good estimate)?
# given the length of 1 array (for each array in a sample), compare it with the reference genome array/ratio: can you say whether or not they're different
    # can we say how many copies they have
# TRID with the suffix _I are labeled indistinguishables; therefore, they may not have accurate read counts/whatever since they may map in a bunch of different places