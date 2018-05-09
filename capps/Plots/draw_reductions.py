import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/reductions.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

pypl.plot(results[0],results[1],label="Tbb parallel_reduce")
pypl.plot(results[0],results[2],label="Homemade reduction")

pypl.xlim(xmin = 0)
pypl.xlabel("Grain size for homemade reduction")
pypl.ylabel("Mean computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Comparison of reductions")
pypl.show()
pypl.savefig("Plots/reductions.jpg");

