import matplotlib.pyplot as plt
from collections import defaultdict
import math
from src.GenomicUtils.LocusFile import LociManager

a=LociManager("/home/avraham/MaruvkaLab/msmutect_runs/data/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed")
sizes = defaultdict(int)
batch_size = 1_000_000
batch=[]
croc=True
while len(batch)==batch_size or croc:
    croc=False
    batch = a.get_batch(batch_size)
    for l in batch:
        if l.locus_length < 131:
            # print(l.chromosome)
            # print(l.start)
            # sizes[l.locus_length]+=1

            sizes[round(l.repeats+0.01)]+=1

print(sizes)
plt.scatter(sizes.keys(), [math.log10(x) for x in sizes.values()])
# plt.xticks(list(sizes.keys()))
plt.yticks(list(range(0, 9)))
plt.ylabel("Number of Loci (Log Base 10 Scale)")
plt.xlabel("Number of repeats (rounded)")
plt.show()
croc=1

