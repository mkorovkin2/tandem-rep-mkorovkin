import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt
from scipy.stats import linregress

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

mlist = list()
slist = list()
tlist = list()
mdict = {}

sim_times = 100
READ_LENGTH = 100
Clist = [70, 71]#355]
Alist = np.arange(READ_LENGTH * 2, READ_LENGTH * 20, READ_LENGTH)

for C in Clist:
    mdict[C] = {}
    for A in Alist:
        mdict[C][A] = (0, 0)

for A in Alist:
    for C in Clist:
        mean, std = sim_adv(A, READ_LENGTH, READ_LENGTH, sim_times, C, 550, 150)
        mdict[C][A] = (mean, std)

fig, axs = plt.subplots(nrows=len(Clist), ncols=1)

print("---")

iter = 0
for C in Clist:
    mAlist = list()
    sAlist = list()
    for A in Alist:
        mAlist.append(mdict[C][A][0])
        sAlist.append(mdict[C][A][1])
    sAlist = np.array(sAlist) * 1.96

    ax = axs[iter]
    ax.errorbar(Alist, mAlist, yerr=[sAlist, sAlist], fmt='o')
    beta1, beta0, rv, pv, stderr = linregress(Alist, mAlist)
    temp_x = np.arange(min(Alist), max(Alist), 10)
    ax.plot(temp_x, [(beta1 * xx + beta0) for xx in temp_x], color='r', alpha=0.7)
    ax.set_title('Coverage is [' + str(C) + "] |  y = " + str(np.round(beta1, decimals=5)) + "x + " + str(np.round(beta0, decimals=5)))
    print(str(np.round(beta1, decimals=7)) + " * A + " + str(np.round(beta0, decimals=7)))
    iter += 1

    #print("\nAdjacency matrix for C=" + str(C))
    #print(find_intersections(mAlist, sAlist, Alist))

plt.ylabel("X/Y ratio")
plt.xlabel("Array length")
#fig.suptitle(str(sim_times) + " simulations for a (READ_LENGTH=" + str(READ_LENGTH) + ")")
plt.show()

print("\nNote on how to interpret the matrices above:\n- each row is a level of A;\n- each column is the corresponding A value;"
      "\n- a 0 indicates that the statistical between A_row and A_column is not significant;"
      "\n- otherwise, it is significant, and the A_column value is listed in the space.")