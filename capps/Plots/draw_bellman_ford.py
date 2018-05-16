import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/bellman_ford.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

pypl.plot(results[0],results[1],label="Doubles")
pypl.plot(results[0],results[2],label="Mpfr")
pypl.plot(results[0],results[3],label="Lazy mpfr, second step of computation")
pypl.plot(results[0],results[4],label="Lazy mpfr, first step with interval arithmetic")
pypl.plot(results[0],results[5],label="Lazy mpfr, total time")

pypl.xlabel("Number of vertices")
pypl.ylabel("Relative mean computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Computation time for Bellman-Ford as a function of the graph size");
pypl.show()
pypl.savefig("Plots/bellmanford.jpg");

