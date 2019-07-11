import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt

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

def sim_adv(A, R, F, sim_times, C, kwargs=""):
    loss = (False, 0, 0)
    gain = (False, 0, 0)
    A_v = (False, None, 0, 0)
    R_v = (False, None, 0, 0)
    F_v = (False, None, 0, 0)
    fragment_imp = (False, 0)
    if len(kwargs) > 0:
        past_kwargs = kwargs
        if "loss" in kwargs:
            to_find = "loss"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            loss = (True, np.float32(followup[0]), np.float32(followup[1]))
        if "gain" in kwargs:
            to_find = "gain"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            gain = (True, np.float32(followup[0]), np.float32(followup[1]))
        if "fragment" in kwargs:
            to_find = "fragment"
            followup = kwargs[(kwargs.index(to_find) + len(to_find) + 1):].strip().split(" ")
            gain = (True, np.int32(followup[0]), np.float32(followup[1]))
        if "variance" in kwargs:
            if "gain" in kwargs:
                to_find = "gain"
                to_find_index = kwargs.index(to_find)
                kwargs = (kwargs[:to_find_index] + kwargs[(to_find_index + len(to_find)):]).strip()
            A_v, R_v, F_v = find_variance(kwargs, A_v, R_v, F_v)

    if A_v[0]:
        A += np.random.normal(0, np.sqrt(A_v[1]), 1)[0]
    if R_v[0]:
        R += np.random.normal(0, np.sqrt(R_v[1]), 1)[0]
    if F_v[0]:
        F += np.random.normal(0, np.sqrt(F_v[1]), 1)[0]

    T = A + 2 * F
    N = np.int32(np.ceil(C * T / R))

    hit_count_list = list()
    out_count_list = list()

    for _ in range(sim_times):
        hit_count = 0
        out_count = 0
        starting_points = list()

        if not fragment_imp[0]:
            for _ in range(N):
                if loss[0]:
                    rand_int = np.random.uniform(0,1)
                    if rand_int < 0.5:
                        randres = np.floor(rand.uniform(0, A + np.int32(0.5 * F))) + loss[1]
                        starting_points.append(randres if randres < A + 0.5 * F else 0)
                    else:
                        randres = np.floor(rand.uniform(0, A + np.int32(0.5 * F))) - loss[1]
                        starting_points.append(randres if randres >= 0 else 0)
                if gain[0]:
                    starting_points.append(np.floor(rand.uniform(gain[1], A + np.int32(0.5 * F)) - gain[1]))
        else:
            random_normal_values = np.random.normal(fragment_imp[1], np.sqrt(fragment_imp[2]), np.int32(N / 2))
            for index in range(N):
                random_index1 = np.floor(rand.uniform(0, A + np.int32(0.5 * F)))
                random_index2 = random_index1 + random_normal_values[index]
                starting_points.append(random_index1)
                starting_points.append(random_index2)

        A_lower = F
        A_upper = (A_lower - 0.25 * F) + A - 0.75 * F - R
        F_lower = 0
        F_lower_bound = F * 0.5
        F_upper = A_upper
        F_upper_bound = F_upper + F * 0.5

        for i in range(N):
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
        ratio_list.append(out_count_list[i] / hit_count_list[i])

    std_return = np.round(np.std(ratio_list), decimals=2)
    mean_return = np.round(np.mean(ratio_list), decimals=2)
    print("A=" + str(A) + ", C=" + str(C) + ", mean=" + str(mean_return) + ", std=" + str(std_return))
    return mean_return, std_return

mlist = list()
slist = list()
tlist = list()

mdict = {}
sim_times = 40

READ_LENGTH = 100
fragment_variance = 50

Clist = [70, 140, 355]
Alist = np.arange(READ_LENGTH * 2, READ_LENGTH * 8)#[400, 550, 700]#[1000, 1500, 4000]#[1000, 1500, 2000]#[500, 750, 1000, 1250, 1500, 1750, 2000]
for C in Clist:
    mdict[C] = {}
    for A in Alist:
        mdict[C][A] = (0, 0)

ordered_mean = list()
ordered_t = list()
test_list = np.arange(READ_LENGTH, 250, fragment_variance)
iter = 0

for t in test_list:
    for A in Alist:
        cmean = 0
        for C in Clist:
            mean, std = sim_adv(A, READ_LENGTH, READ_LENGTH, sim_times, C, kwargs="fragment " + str(t) + " " + str(fragment_variance))
            mdict[C][A] = (mean, std)
            cmean += mean
        cmean = np.round(cmean / len(Clist), decimals=2)
        ordered_mean.append(cmean)
        ordered_t.append(A)
        print("> Progress: " + str(100 * np.round(iter / (len(test_list) * len(Alist)), decimals=2)) + "%")
        iter += 1

for index in range(len(Alist)):
    plt.scatter(test_list, ordered_mean[index::len(Alist)])

plt.title("READ_LENGTH=" + str(READ_LENGTH) + " | " + str(sim_times) + " iterations per simulation")
plt.ylabel("X/Y ratio")
plt.xlabel("t value (the intermediate fragment length | std.dev=" + str(np.round(np.sqrt(fragment_variance), decimals=2)) + ")")
plt.show()

# Fragment size from left end to right end; read lengths [100, 150, 250]; separation of coverages; increase number of simulations to 100;
#       allow fragment length to be smaller than 2R; fragment lengths 400-700 (550~150);
# Graph of paired reads with error bars; -> add regression line; regression on 95% confidence intervals
# Describe simulations in text on a document
# loss: losing data from re-mapping reads (chance of loss?)

# ultimately: how different do ratios have to be so that you can tell that array length changed; can you tell the length of the array (good estimate)?
# given the length of 1 array (for each array in a sample), compare it with the reference genome array/ratio: can you say whether or not they're different
    # can we say how many copies they have
# TRID with the suffix _I are labeled indistinguishables; therefore, they may not have accurate read counts/whatever since they may map in a bunch of different places