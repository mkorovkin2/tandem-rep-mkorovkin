import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd

def sim_adv(A, R, F, sim_times, C, frag_mean, frag_std):
    T = A + 2 * F * 0.75
    N = np.int32(np.ceil(C * T / R))

    # Lists used for tracking read hits/misses
    hit_count_list = list()
    out_count_list = list()

    # Simulate sim_times number of times
    for _ in range(sim_times):
        hit_count = 0
        out_count = 0

        # List starting points for fragments
        starting_points = list()

        # Generate a random set of simulated DNA read fragments
        random_normal_values = np.random.normal(frag_mean, frag_std, np.int32(N / 2))

        # Iterate through fragment array, adding 2 reads for each fragment to STARTING_POINTS
        for index in range(np.int32(N / 2)):
            random_index1 = np.floor(rand.uniform(0, A + np.int32(0.5 * F)))
            random_index2 = random_index1 + random_normal_values[index]

            starting_points.append(random_index1)
            starting_points.append(random_index2)

        # Calculate acceptable boundaries for hits/misses
        A_lower = F * 0.75 # A_lower = F
        A_upper = (A_lower - 0.25 * F) + A - R # A_upper = (A_lower - 0.25 * F) + A - 0.75 * F - R
        F_lower = 0
        F_lower_bound = F * 0.5
        F_upper = A_upper
        F_upper_bound = F_upper + F * 0.5

        # Iterate through starting points; identify hits and misses
        for index in range(np.int32(N / 2)):
            x = starting_points[index]
            if (F_lower <= x and F_lower_bound > x) or (F_upper < x and F_upper_bound >= x):
                hit_count += 1
            elif (A_lower <= x and A_upper >= x):
                out_count += 1

        # Append total hit count
        hit_count_list.append(hit_count)
        # Append total miss count
        out_count_list.append(out_count)

    # Total humbers of hits/misses
    sum_hit = 0
    sum_out = 0

    # List of ratios (hit/miss)
    ratio_list = list()

    for i in range(sim_times):
        if hit_count_list[i] > 0:
            # Sum the number of hits and misses
            sum_hit += hit_count_list[i]
            sum_out += out_count_list[i]
            # Append the hit/miss ratio to the RATIO_LIST
            ratio_list.append(out_count_list[i] / hit_count_list[i])

    # Calculate mean/standard deviation of RATIO_LIST
    std_return = np.round(np.std(ratio_list), decimals=2)
    mean_return = np.round(np.mean(ratio_list), decimals=2)

    # Print out statistics from the simulations
    print("A=" + str(A) + ", C=" + str(C) + ", mean=" + str(mean_return) + ", std=" + str(std_return))

    # Return all useful information
    return mean_return, std_return

def sim_adv_new(read_length, A, F, frag_mean, frag_std, sim_times):
    # Values to calculate
    C = 355#100
    g = 10000#read_length * 51
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

        # List starting points for fragments
        starting_points = list()

        # Generate a random set of simulated DNA read fragments
        random_normal_values = np.random.normal(frag_mean, frag_std, number_fragments_from_g)
        random_starting_points = np.random.uniform(0, (G - frag_mean * 3), number_fragments_from_g)

        #random_array_location = np.random.uniform(np.int32(A / 4), np.int32(G - (A + (A / 4))), 1)[0]

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

def init_sim(R):
    Alist = np.arange(R * 3, R * 40, R)
    mAlist = list()
    sAlist = list()
    for A in Alist:
        mean, std = sim_adv_new(R, A, np.int32(R * 0.5), 550, 150, 5000)
        mAlist.append(mean)
        sAlist.append(std * 1.96)
        print(R, A)

    plt.errorbar(Alist, mAlist, yerr=[sAlist, sAlist])
    print(linregress(Alist, mAlist)[0:2])
    plt.title("READ_LENGTH=" + str(R) + " | Graph of array length vs. X/Y ratio")
    plt.ylabel("X/Y ratio")
    plt.xlabel("Array length")
    plt.show()

    return np.array(mAlist), np.array(sAlist)

READ_LENGTH = [100, 148, 250]

mA0, sA0 = init_sim(READ_LENGTH[0])
mA1, sA1 = init_sim(READ_LENGTH[1])
mA2, sA2 = init_sim(READ_LENGTH[2])

overall_df = pd.DataFrame(data={"mA100": mA0,
                                "sA100": sA0,
                                "mA148": mA1,
                                "sA148": sA1,
                                "mA250": mA2,
                                "sA250": sA2
                                }
                          )
overall_df.to_csv("/Users/mkorovkin/Desktop/marzd/output40_coverage355_00.csv")