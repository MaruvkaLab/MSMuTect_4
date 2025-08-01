from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import math, os, glob, shutil
import numpy as np

def dominant_percentage(s):
    a={"A": 0, "C": 0, "G": 0, "T": 0}
    for c in s:
        a[c]+=1
    return max(a.values())/len(s)

def analyze_file(fp: str):
    df = pd.read_csv(fp, delimiter="\t")
    ret = np.zeros((15, 20))
    for index, row in df.iterrows():
        # print(row.NORMAL_MOTIF_REPEATS_1)
        # print(row.TUMOR_SUPPORTING_READS_2)
        if row.CALL=='M':

            dom = dominant_percentage(row.PATTERN)
            ret[len(row.PATTERN) - 1, int((round(dom, 2) - 0.05) * 20)] += 1
    plt.imshow(ret)
    plt.show()
    return ret

def analyze_directory(directory: str, out_file_name):
    all_files = glob.glob(os.path.join(directory, "*.tsv.gz"))
    out_npy = np.zeros((len(all_files), 15, 20))
    for i, f in enumerate(all_files):
        os.system(f"gunzip -k {f}")
        tmp_file = f[:-3]
        current_dist = analyze_file(tmp_file)
        out_npy[i] = current_dist
        os.remove(tmp_file)
    np.save(out_file_name, out_npy)


# analyze_directory("/home/avraham/MaruvkaLab/Texas/entropy_analysis/mss", "mss")
analyze_directory("/home/avraham/MaruvkaLab/Texas/entropy_analysis/msi", "msi")