import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/lazympsseqalt.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

#pypl.plot(results[0],results[1],label="Doubles")
#pypl.plot(results[0],results[2],label="Summation with superaccumulators")
#pypl.plot(results[0],results[3],label="Mps with superaccumulators")
#pypl.plot(results[0],results[4],label="Mps with mpfr")
pypl.plot(results[0],results[5],label="Interval over-approximation")
pypl.plot(results[0],results[6],label="Reverse mps processing")
pypl.plot(results[0],results[7],label="Reverse pos processing")
pypl.plot(results[0],results[8],label="Exact mps with superacc")
pypl.plot(results[0],results[9],label="Exact pos with superacc")
#pypl.plot(results[0],results[10],label="Exact mps with mpfr")
#pypl.plot(results[0],results[11],label="Exact pos with mpfr")
pypl.plot(results[0],results[12],label="Lazy mps with superacc")
pypl.plot(results[0],results[13],label="Lazy pos with superacc")
#pypl.plot(results[0],results[14],label="Lazy mps with mpfr")
#pypl.plot(results[0],results[15],label="Lazy pos with mpfr")


pypl.xlabel("Dynamic range")
pypl.ylabel("Relative mean computation time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Computation time for mps as a function of the dynamic range");
pypl.show()
pypl.savefig("Plots/lazympsseq.jpg");

