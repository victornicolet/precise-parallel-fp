import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/parallel_mss.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

pypl.plot(results[0],results[1],label="Sequential doubles")
pypl.plot(results[0],results[2],label="Parallel doubles")
pypl.plot(results[0],results[3],label="Dynamic lazy computation")
#pypl.plot(results[0],results[4],label="Lazy computation with superaccumulators")
#pypl.plot(results[0],results[5],label="Lazy computation with superaccumulators, optional optimization")

pypl.xlabel("Dynamic range")
pypl.ylabel("Relative mean computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Computation time for parallel mss as a function of the dynamic range");
pypl.show()
pypl.savefig("Plots/lazymsspar.jpg");

