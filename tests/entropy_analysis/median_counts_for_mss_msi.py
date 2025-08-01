import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
# a=pd.read_csv("/home/avraham/MaruvkaLab/p/n/mss_counts.tsv", delimeter="\t")


import pandas as pd

df = pd.read_csv("/home/avraham/MaruvkaLab/p/n/msi_counts.tsv", sep="\t")
median_values_msi = df.iloc[:, :16].median().to_dict()
df = pd.read_csv("/home/avraham/MaruvkaLab/p/n/mss_counts.tsv", sep="\t")
median_values_mss = df.iloc[:, :16].median().to_dict()

num_muts_mss = np.array([median_values_mss[f"Motif_size{i}_mut_count"] for i in range(1, 16)])

num_muts_msi = np.array([median_values_msi[f"Motif_size{i}_mut_count"] for i in range(1, 16)])

plt.yscale('log')  # This is the key line

num_loci = np.array([1057987, 1706448, 21849768, 2282940, 417426, 193841, 5156, 68481, 9959, 20236, 40467, 8162, 17095, 5339, 5853])
plt.plot(np.arange(1,16), num_muts_msi/num_loci, color="green", label="MSI")
plt.plot(np.arange(1,16), num_muts_mss/num_loci,  color="red", label="MSS")
plt.legend()
plt.ylabel("Median Mutation Frequency")
plt.xlabel("Motif Length")
# plt.yticks(np.arange(0, 0.25, 0.02))
plt.show()


plt.plot(np.arange(1,16), np.log10(num_muts_msi), color="green", label="MSI")
plt.plot(np.arange(1,16), np.log10(num_muts_mss),  color="red", label="MSS")
plt.show()


# Data
x_values = [8, 7]
heights = [0.99, 0.98]
colors = ['blue', 'red']

# Create bar graph
plt.bar(x_values, heights, color=colors)
plt.xticks(list(range(5, 10)))

# Optional: label axes and show the plot
plt.xlabel("Repeats")
plt.ylabel("Allele Fraction")
plt.legend()
plt.title("Reversion to Reference Filtered")
legend_elements = [
    Patch(facecolor='blue', label='Normal'),
    Patch(facecolor='red', label='Tumor (and Reference)')
]
plt.legend(handles=legend_elements, title='Legend')
plt.show()



#plt.hist([8 for i in range(15)]+[6 for j in range(12)], width=0.9, align='left')
x={6:0.45, 8:0.55}
keys = list(x.keys())
values = list(x.values())
plt.hist(keys, values, width=0.8, align='left')
#
# plt.xticks(list(range(5, 15)))
plt.xlabel("Number of Repeats")
plt.ylabel("Allele Frequency")
plt.show()


# data = {8: 0.46, 6: 0.38, 5:0.15}
data = {8: 0.99}

keys = list(data.keys())
values = list(data.values())

plt.bar(keys, values, width=0.8, align='center')  # you can adjust width to make bars "fatter"
plt.xlabel("Number of Repeats")
plt.ylabel("Allele Frequency")
# plt.title('Custom Value Histogram')
plt.xticks(keys)  # ensures both 7 and 8 show up cleanly
plt.show()