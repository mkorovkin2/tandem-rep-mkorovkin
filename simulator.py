import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt

#A = 1000
#R = 250
#F = 125
#C = 80

def find_intersections(mAlist, sAlist, Alist):
    adjacency_matrix = np.zeros(shape=(len(mAlist), len(mAlist)))
    for mi in range(len(mAlist)): # row: mi, column: smi (should be symmetric)
        std = sAlist[mi]
        mean = mAlist[mi]
        for smi in range(len(mAlist)):
            substd = sAlist[smi]
            submean = mAlist[smi]
            if (mean + std < submean - substd) or (mean - std > submean + substd):
                adjacency_matrix[mi][smi] = Alist[smi]
                #print("testing", Alist[mi], "against", Alist[smi], "result:", Alist[mi])
            else:
                adjacency_matrix[mi][smi] = 0
                #print("testing", Alist[mi], "against", Alist[smi])
    return adjacency_matrix

def sim(A, R, F, sim_times, C=80):
    T = A + 2 * F
    N = np.int32(np.ceil(C * T / R))

    hit_count_list = list()
    out_count_list = list()

    # Simulate 100 times
    for _ in range(sim_times):
        hit_count = 0
        out_count = 0
        starting_points = list()

        for _ in range(N):
            starting_points.append(np.floor(rand.uniform(0, A)))

        # define boundaries
        A_lower = F
        A_upper = (A_lower - 0.25 * F) + A - 0.75 * F - R
        F_lower = 0
        F_lower_bound = F * 0.5
        F_upper = A_upper
        F_upper_bound = F_upper + F * 0.5



        for i in range(N):
            x = starting_points[i]
            if (F_lower <= x and F_lower_bound > x) or (F_upper < x and F_upper_bound >= x):#(x < Fs) or ((x > (A + F) - R) and (x < A + Fs)):
                hit_count += 1
            elif (A_lower <= x and A_upper >= x):#(x > Fin) and (x < (A + Fin) - R):
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

    #print("Mean ratio:", np.mean(ratio_list))
    #print("Median ratio:", np.median(ratio_list))

    plt.hist(ratio_list, bins=20)
    plt.title("A=" + str(A) + ", R=" + str(R) + " F=" + str(F) + " C=" + str(C) + " | Mean ratio: " + str(np.round(np.mean(ratio_list), decimals=2)) + " | Std.Dev: " + str(np.round(np.std(ratio_list), decimals=2)))
    plt.xlabel("Ratio of internal to hit (internal / hit)")
    plt.ylabel("Frequency of ratio")
    #plt.savefig("/Users/mkorovkin/Desktop/TRImages/figure" + str(id) + ".png")
    plt.clf()

    std_return = np.round(np.std(ratio_list), decimals=2)
    mean_return = np.round(np.mean(ratio_list), decimals=2)
    print("A=" + str(A) + ", C=" + str(C) + ", mean=" + str(mean_return) + ", std=" + str(std_return))
    return mean_return, std_return

#sim(1000, 250, 125)

mlist = list()
slist = list()
tlist = list()

mdict = {}
sim_times = 100

Clist = [70, 140, 355] # [70, 140, 355]
Alist = [500, 750, 1000, 1250, 1500, 1750, 2000] # [500, 750, 1000, 1250, 1500, 1750, 2000]:
for C in Clist:
    mdict[C] = {}
    for A in Alist:
        mdict[C][A] = (0, 0)

for A in Alist:
    for C in Clist:
        mean, std = sim(A, 100, 100, sim_times, C)#sim(A, 250, 200, sim_times, C) # sim(A, 100, 50, id, C)
        mdict[C][A] = (mean, std)
        #mlist.append(mean)
        #slist.append(std)
        #tlist.append((A, C))

# Now switch to a more OO interface to exercise more features.
fig, axs = plt.subplots(nrows=len(Clist), ncols=1)
print(axs)

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
    ax.set_title('Graph C=' + str(C))
    iter += 1

    print("Adjacency matrix for C=", str(C))
    print(find_intersections(mAlist, sAlist, Alist))
    print("Note: how to interpret: each row is a level of A; each column is the corresponding A value; a 0 indicates that the statistical between A_row and A_column is not significant; otherwise, it is significant, and the A_column value is listed in the space.")

plt.ylabel("X/Y ratio")
plt.xlabel("A value")
plt.show()

# Summary: for every 250 increase in A, mean of ratio increases by 1. Changes in C have no effect on the ratio.