from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import math, os, glob, shutil, pickle
import numpy as np


def analyze_file(fp: str):
    df = pd.read_csv(fp, delimiter="\t")
    ret = defaultdict(int)
    for index, row in df.iterrows():
        if row.CALL=='M':
            ret[row.PATTERN]+=1
    return ret

def analyze_directory(directory: str, out_file_dir: str):
    all_files = glob.glob(os.path.join(directory, "*.tsv.gz"))
    for i, f in enumerate(all_files):
        os.system(f"gunzip -k {f}")
        tmp_file = f[:-3]
        case_id = os.path.basename(tmp_file)[:-20]

        current_pattern_dict = analyze_file(tmp_file)
        os.remove(tmp_file)
        with open(f'{os.path.join(out_file_dir, case_id)}.pkl', 'wb') as f:
            pickle.dump(current_pattern_dict, f)


analyze_directory("/home/avraham/MaruvkaLab/Texas/entropy_analysis/mss", "/home/avraham/MaruvkaLab/Texas/entropy_analysis/mss_dict")
analyze_directory("/home/avraham/MaruvkaLab/Texas/entropy_analysis/msi", "/home/avraham/MaruvkaLab/Texas/entropy_analysis/msi_dict")