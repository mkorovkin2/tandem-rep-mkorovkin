import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd

def sim_adv_new(read_length, A, F, frag_mean, frag_std, sim_times, gain_loss):
    # Values to calculate
    C = 355
    g = 10000
    G = g + A + g
    number_fragments = np.int32(np.ceil(C * (A) / frag_mean))
    number_fragments_from_g = np.int32(C * G / frag_mean)

    # Lists used for tracking read hits/misses
    hit_count_list = list()
    out_count_list = list()

    random_array_location = g

    # Simulate sim_times number of times
    for _ in range(sim_times):
        hit_count = 0
        out_count = 0

        # Generate a random set of simulated DNA read fragments
        random_normal_values = np.random.normal(frag_mean, frag_std, number_fragments_from_g)
        random_starting_points = np.random.uniform(0, (G - frag_mean * 3), number_fragments_from_g)

        A_lower = F * 0.75 + random_array_location + (F * 0.25)
        A_upper = A_lower + A - read_length - (F * 0.25)
        F_lower = 0 + random_array_location
        F_lower_bound = F * 0.5 + random_array_location
        F_upper = A_upper + F * 0.25
        F_upper_bound = A_upper + read_length# F_upper + F * 0.75

        # print(F_lower, F_lower_bound, A_lower, A_upper, F_upper, F_upper_bound)

        # Iterate through fragment array, adding 2 reads for each fragment to STARTING_POINTS
        for index in range(len(random_normal_values)):
            x1 = random_starting_points[index]
            x2 = random_starting_points[index] + random_normal_values[index] - read_length

            for x in [x1, x2]:
                if (F_lower <= x and F_lower_bound > x) or (F_upper < x and F_upper_bound >= x):
                    hit_count += 1
                elif (A_lower <= x and A_upper >= x):
                    out_count += 1

        if hit_count > 0:
            hit_count_list.append(hit_count)
            out_count_list.append(out_count)

    # Return all useful information
    hit_over_out = np.array(out_count_list) / np.array(hit_count_list)
    return hit_over_out.mean(), hit_over_out.std()

def init_sim(R, gain_loss):
    Alist = np.arange(R * 3, R * 40, R)
    mAlist = list()
    sAlist = list()
    for A in Alist:
        mean, std = sim_adv_new(R, A, np.int32(R * 0.5), 550, 150, 100, gain_loss)
        mAlist.append(mean)
        sAlist.append(std * 1.96)
        #print(R, A)

    plt.errorbar(Alist, mAlist, yerr=[sAlist, sAlist])
    print(linregress(Alist, mAlist)[0:2])
    plt.title("READ_LENGTH=" + str(R) + " | Graph of array length vs. X/Y ratio")
    plt.ylabel("X/Y ratio")
    plt.xlabel("Array length")
    plt.show()

    C = 3.8
    print("eq00 =", linregress(Alist, mAlist)[0:2])
    print("eq01 =", linregress(Alist, np.array(mAlist) * (C + 0.5) / C)[0:2])
    print("eq0n1 =", linregress(Alist, np.array(mAlist) * (C - 0.5) / C)[0:2])
    print("eq11 =", linregress(Alist, np.array(mAlist) * (C + 1) / C)[0:2])
    print("eqn11 =", linregress(Alist, np.array(mAlist) * (C - 1) / C)[0:2])

    print("std =", linregress(Alist, sAlist)[0:2])

    return np.array(mAlist), np.array(sAlist)

READ_LENGTH = [100, 148, 250]
gain_loss = 0

mA0, sA0 = init_sim(READ_LENGTH[0], gain_loss)
mA1, sA1 = init_sim(READ_LENGTH[1], gain_loss)
mA2, sA2 = init_sim(READ_LENGTH[2], gain_loss)

exit(0)

overall_df = pd.DataFrame(data={"mA100": mA0,
                                "sA100": sA0,
                                "mA148": mA1,
                                "sA148": sA1,
                                "mA250": mA2,
                                "sA250": sA2
                                }
                          )

overall_df.to_csv("/Users/mkorovkin/Desktop/marzd/output_statistics_"
                  + str(np.int32(np.ceil(gain_loss))) + str(np.int32(np.ceil(gain_loss)))
                  + ".csv")