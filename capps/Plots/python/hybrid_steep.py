import csv
import matplotlib.pyplot as pypl
import numpy as np

results = []
with open("Plots/csv/steep_hybrid.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        try:
            results.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
curves = np.array(results).transpose()
print curves

results2 = []
with open("Plots/csv/mss_hybrid_final.csv", "rb") as csvfile2:
    resultsreader2 = csv.reader(csvfile2)
    for row in resultsreader2:
        try:
            results2.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
curves2 = np.array(results2).transpose()
print curves2


pypl.figure(1)
pypl.subplot(221)
pypl.plot(curves2[0],curves2[4]/curves2[1],label="Double Parallel Reduce")
pypl.plot(curves2[0],curves2[4]/curves2[2],label="Interval Arithmetic Filtering")
pypl.plot(curves2[0],curves2[4]/curves2[4],label="Sequential Implementation")
pypl.fill_between(curves2[0],curves2[4]/curves2[2],curves2[4]/curves2[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves2[0],curves2[4]/curves2[3],curves2[4]/curves2[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')

pypl.xlim(3*10**5,10**9)
pypl.ylim(0.5,5)
pypl.xlabel("Array Size")
pypl.ylabel("Relative Throughput")
pypl.xscale('log')
pypl.yticks(np.arange(0.5,8.5,1))
pypl.title('Mss, machine A')
#pypl.legend(loc=2, prop={'size': 5})


pypl.subplot(222)
pypl.plot(curves[0],curves[4]/curves[1],label="Double Parallel Reduce")
pypl.plot(curves[0],curves[4]/curves[2],label="Interval Arithmetic Filtering")
pypl.plot(curves[0],curves[4]/curves[4],label="Sequential Implementation")
pypl.fill_between(curves[0],curves[4]/curves[2],curves[4]/curves[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves[0],curves[4]/curves[3],curves[4]/curves[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')

pypl.xlim(3*10**5,10**9)
pypl.ylim(0.5,4)
pypl.xlabel("Array Size")
#pypl.ylabel("Relative Throughput")
pypl.xscale('log')
pypl.yticks(np.arange(0.5,4,0.5))
pypl.title('Steep, machine A')

pypl.legend(loc=9, bbox_to_anchor = (0.5,-0.1), prop={'size': 10})

pypl.suptitle("Relative Throughput as a Function of Input Size");
#pypl.show()
pypl.savefig("Plots/figures/hybridsteepfinal.jpg");
