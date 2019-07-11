import numpy as np

content = []
with open("/Users/mkorovkin/Desktop/bp250_mk2.o2625098") as f:#bp100.coverage.bed")as f:
    #
    #
    #
    # WORKS WITH THE .BED FILE BUT DOESN'T WORK WITH THIS - FLANKS ARE OFF?
    #
    #
    for line in f:
        content.append(line.strip().split())

content2 = []
with open("/Users/mkorovkin/Desktop/bp150_mk2.o2625011") as f:
    for line in f:
        content2.append(line.strip().split())

content3 = []
with open("/Users/mkorovkin/Desktop/bp100_mk2.o2624915") as f:
    for line in f:
        content3.append(line.strip().split())
#bp100_mk2.o2624915
#bp150_mk2.o2625011
#bp250_mk2.o2625098

def analyze(content):
    flank1_read = {}
    flank1_bp = {}
    flank2_read = {}
    flank2_bp = {}
    array_read = {}
    array_bp = {}
    unique_keys = list() ## readcount coveredbp arraysize coverage .>
    remove_keys = list()

    for item in content:
        try:
            sub = item[3]
            appended = False
            print(item)
            key = str(sub.split("_")[0])
            if len(str(sub)) > 10 and ("flank1" in sub):
                flank1_read[key] = int(item[4])
                flank1_bp[key] = int(item[5])
            elif len(str(sub)) > 10 and ("flank2" in sub):
                flank2_read[key] = int(item[4])
                flank2_bp[key] = int(item[5])
            else:
                array_read[key] = int(item[4])
                array_bp[key] = int(item[5])
            unique_keys.append(key)
        except:
            print("Error")

    print(len(unique_keys))
    sum_ratio = 0
    count = 0
    ratio_dict = {}
    dlist = list()
    upper = 1500
    lower = 0
    for key in unique_keys:
        #print(flank1_read[key], flank2_read[key], array_read[key])
        if (flank1_read[key] + flank2_read[key]) > 0 and array_read[key] < upper and array_read[key] > lower:
            ratio = array_read[key] / (flank1_read[key] + flank2_read[key])
            if ratio > 10:
                pass#print(array_read[key], flank1_read[key], flank2_read[key])
            dlist.append(ratio)
            ratio_dict[ratio] = (array_bp[key] + flank1_bp[key] + flank2_bp[key])
            sum_ratio += ratio
            count += 1
    return ratio_dict

ratio_dict1 = analyze(content)
ratio_dict2 = analyze(content2)
ratio_dict3 = analyze(content3)

import matplotlib.pyplot as plt
#plt.hist(dlist, bins=50, range=[0, 15])
plt.scatter(list(ratio_dict1.values()), list(ratio_dict1.keys()))               #250
plt.scatter(list(ratio_dict2.values()), list(ratio_dict2.keys()), alpha=0.7)    #150
plt.scatter(list(ratio_dict3.values()), list(ratio_dict3.keys()), alpha=0.4)    #100
plt.xlabel("Total length of analyzed sequences")
plt.ylabel("Ratio X/Y")
plt.show()

#from scipy.stats import linregress
#print(linregress(list(ratio_dict.values()), list(ratio_dict.keys())))

# With smaller read lengths (<1000), the ratio seems to be closer to 1.5ish;
# Medium size read lengths (<2500), the ratio tends to be 3.5ish;
# The whole set ratio is around 10.4


# QUESTIONS:
# 0. How do I know which read lengths are on each of these files?
# 1. What is the "simulation" that you sent (i.e. bed file)?
# 2. How do I actually log into the database? Which files should I use/take?
# 3. What exactly is it that I'm trying to do? I'm a little confused on this...