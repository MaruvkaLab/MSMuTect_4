from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

def dominant_percentage(s):
    a={"A": 0, "C": 0, "G": 0, "T": 0}
    for c in s:
        a[c]+=1
    return max(a.values())/len(s)


hetero = np.load("../positive_and_negative_model/dom_percentages_loci.npy")
plt.imshow(np.log(hetero))
plt.show()
fp = "/home/avraham/MaruvkaLab/msmutect_runs/negatives/m.tsv"
df = pd.read_csv(fp, sep='\t')
croc=1

sl = [[] for i in range(20)]
x = defaultdict(int)
ret = np.zeros((15, 20))

for index, row in df.iterrows():
    # print(row.NORMAL_MOTIF_REPEATS_1)
    # print(row.TUMOR_SUPPORTING_READS_2)
    dom = dominant_percentage(row.PATTERN)
    ret[len(row.PATTERN) - 1, int((round(dom, 2) - 0.05) * 20)] += 1

    # rs2 = row.TUMOR_SUPPORTING_READS_2
    # if math.isnan(rs2):
    #     continue
    # rs2 = int(rs2)
    # if row.NORMAL_MOTIF_REPEATS_1==row.TUMOR_MOTIF_REPEATS_1 and rs2>=5 and math.isnan(row.NORMAL_MOTIF_REPEATS_2):
    #     # print(rs2)
    #     if False and row.TUMOR_FRACTION_2<0.13:
    #         print(f"N1: {row.NORMAL_MOTIF_REPEATS_1}")
    #         print(f"N2: {row.NORMAL_MOTIF_REPEATS_2}")
    #         print(f"T1: {row.TUMOR_MOTIF_REPEATS_1}")
    #         print(f"T2: {row.TUMOR_MOTIF_REPEATS_2}")
    #         print(f"NR: {row.NORMAL_SUPPORTING_READS_1}")
    #         print(f"TR: {row.TUMOR_SUPPORTING_READS_1}")
    #         print(f"AIC TUMOR NORMAL: {row.AIC_TUMOR_NORMAL-row.AIC_TUMOR_TUMOR}")
    #         print(f"AIC NORMAL TUMOR: {row.AIC_NORMAL_TUMOR-row.AIC_NORMAL_NORMAL}")
    #
    #         # print(row.NORMAL_SUPPORTING_READS_1)
    #         # print(row.NORMAL_FRACTION_1)
    #         # print(row.TUMOR_FRACTION_1)
    #         # print(row.TUMOR_FRACTION_2)
    #         croc=1
    #
    #     x[round(row.TUMOR_FRACTION_2, 2)] += 1
    #     sl[int(round(row.TUMOR_FRACTION_2, 2)*20)].append(row.NORMAL_SUPPORTING_READS_1)
croc=1
plt.scatter(x.keys(), x.values())
plt.ylabel("Number of mutations")
plt.xlabel("Allele Fraction of tumor mutation")
plt.show()

plt.scatter(np.arange(0, 1, 0.05), [np.median(np.array(c)) for c in sl])
plt.ylabel("Number of reads")
plt.xlabel("Allele Fraction of tumor mutation")
plt.show()
