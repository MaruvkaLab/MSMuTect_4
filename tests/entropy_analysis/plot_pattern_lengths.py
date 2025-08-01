from typing import List, Set
import matplotlib.pyplot as plt
import math, os, glob, shutil, pickle

from jupyter_core.version import pattern

FREQUENCY_THRESHOLD = 1e-9

def load_dict(dict_path):
    with open(dict_path, 'rb') as f:
        locus_counts = pickle.load(f)
    return locus_counts


def get_patterns(filtered_dicts) -> Set[str]:
    ret = set()
    for f in filtered_dicts:
        for key, val in f.items():
            ret.add(key)
    return ret

def safe_retrieve(d, pattern):
    if pattern in d:
        return d[pattern]
    else:
        return 0

def plot_dicts(filtered_dicts, title):
    patterns = get_patterns(filtered_dicts)
    x_indices = range(len(patterns))

    # Plot each dictionary as a separate scatter series
    for d in filtered_dicts:
        y_values = [safe_retrieve(d, base) for base in patterns]
        plt.scatter(x_indices, y_values)

    # Labeling the x-axis with base letters
    plt.xticks(x_indices, patterns)
    plt.xlabel("Base")
    plt.ylabel("Count")
    plt.title(title)
    # plt.savefig(title, dpi=300)  # You can change the filename and format

    plt.show()

def plot_pattern_length(pattern_length, dicts_dir, title):
    dicts = glob.glob(os.path.join(dicts_dir, "*.pkl"))
    locus_counts = load_dict('../positive_and_negative_model/locus_counts.pkl')
    loaded_dicts = [load_dict(d) for d in dicts]
    filtered_dicts = [{k: v/locus_counts[k] for k, v in d.items() if len(k) == pattern_length and v/locus_counts[k] > FREQUENCY_THRESHOLD} for d in loaded_dicts]
    plot_dicts(filtered_dicts, title+f": {pattern_length}")



def plot_dict_scatter_plot():
    for pattern_length in range(1, 9):
        plot_pattern_length(pattern_length, "/home/avraham/MaruvkaLab/Texas/entropy_analysis/msi_dict", "msi")
        # plot_pattern_length(pattern_length, "/home/avraham/MaruvkaLab/Texas/entropy_analysis/mss_dict", "mss")

if __name__ == '__main__':
    plot_dict_scatter_plot()