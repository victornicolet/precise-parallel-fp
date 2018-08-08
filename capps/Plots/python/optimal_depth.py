import csv
import matplotlib.pyplot as pypl
import numpy as np


results = []
with open("Plots/csv/optimal_depth.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        try:
            results.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
curves = (np.array(results[1:])).transpose()
print curves

pypl.plot(curves[0],curves[1])

pypl.xlabel("Array Size")
pypl.ylabel("Optimal Depth Threshold")

pypl.xscale('log')
pypl.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
pypl.xlim(3*10**3,3*10**8)

pypl.title("Optimal Depth Threshold as a Function of Array Size")
#pypl.show()
pypl.savefig("Plots/figures/optimalDepth.jpg");
