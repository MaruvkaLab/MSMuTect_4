import matplotlib.pyplot as plt
from collections import defaultdict
import math
from src.GenomicUtils.LocusFile import LociManager
import numpy as np

def dominant_percentage(s):
    a={"A": 0, "C": 0, "G": 0, "T": 0}
    for c in s:
        a[c]+=1
    return max(a.values())/len(s)
a=LociManager("/home/avraham/MaruvkaLab/msmutect_runs/data/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed")
# sizes = defaultdict(int)
# num_loci = defaultdict(int)
sizes = [[] for i in range(1, 16)]
batch_size = 1_000_000
batch=[]
croc=True
x=defaultdict(int)
ret = np.zeros((15, 20))
while len(batch)==batch_size or croc:
    croc=False
    batch = a.get_batch(batch_size)
    for l in batch:
        if l.locus_length < 131:
            # print(l.chromosome)
            # print(l.start)
            # sizes[l.locus_length]+=1
            # if l.repeat_length > 15:
            #     print(l.repeat_length)
            # sizes[l.repeat_length-1].append(l.locus_length)
            # dom = dominant_percentage(l.pattern)
            # ret[len(l.pattern)-1, int((round(dom, 2)-0.05) * 20)]+=1
            x[l.pattern]+=1
print(x)
croc=1
np.save("../positive_and_negative_model/dom_percentages_loci.npy", ret)
# sizes_np = [np.array(c) for c in sizes]
# print([np.median(c) for c in sizes_np])
# plt.scatter(sizes.keys(), [math.log10(x) for x in sizes.values()])
# # plt.xticks(list(sizes.keys()))
# plt.yticks(list(range(0, 9)))
# plt.ylabel("Number of Loci (Log Base 10 Scale)")
# plt.xlabel("Number of repeats (rounded)")
# plt.show()
# croc=1

