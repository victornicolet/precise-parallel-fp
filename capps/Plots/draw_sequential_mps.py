import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/plot2.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

pypl.plot(results[0],results[1],label="Doubles")
pypl.plot(results[0],results[2],label="Summation with superaccumulators")
pypl.plot(results[0],results[3],label="Mps with superaccumulators")
pypl.plot(results[0],results[4],label="Lazy computation with superaccumulators")
pypl.plot(results[0],results[5],label="Lazy computation with superaccumulators, optional optimization")

pypl.xlabel("Dynamic range")
pypl.ylabel("Relative mean computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Computation time for mps as a function of the dynamic range");
pypl.show()
pypl.savefig("Plots/lazympsseq.jpg");

