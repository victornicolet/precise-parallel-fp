import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/csv/mss_hybrid.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        try:
            results.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
ref = results[0][0]
curves = (np.array(results[1:])).transpose()
print ref
print curves

pypl.plot(curves[0],curves[1]/ref,label="Lazy computation")
pypl.plot(curves[0],curves[2]/ref,label="Hybrid reduction")

pypl.xlabel("Static depth")
pypl.ylabel("Computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Computation time for parallel mss as a function of the static depth");
pypl.show()
pypl.savefig("Plots/figures/hybridmss.jpg");
