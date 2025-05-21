from scipy.stats import poisson
import matplotlib.pyplot as plt

plots=[]
depth=60
read_length=151
for locus_length in range(0, 130):
    possible_bases = (read_length-locus_length-20)
    # print(possible_bases)
    poisson_lambda = (depth/read_length)*possible_bases
    print(poisson_lambda)
    cdf = [poisson.cdf(x, poisson_lambda) for x in range(200)]
    plots.append(cdf)

for p in plots:
    plt.plot(p)
plt.show()