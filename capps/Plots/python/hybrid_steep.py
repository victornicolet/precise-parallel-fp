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
with open("Plots/csv/mss_hybrid_final_with_recomp.csv", "rb") as csvfile2:
    resultsreader2 = csv.reader(csvfile2)
    for row in resultsreader2:
        try:
            results2.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
curves2 = np.array(results2).transpose()
print curves2

results3 = []
with open("Plots/csv/steep16.csv", "rb") as csvfile2:
    resultsreader2 = csv.reader(csvfile2)
    for row in resultsreader2:
        try:
            results3.append(map(float, row))
        except ValueError,e:
            print row


# Extract result
curves3 = np.array(results3).transpose()
print curves3

results4 = []
with open("Plots/csv/mss16.csv", "rb") as csvfile2:
    resultsreader2 = csv.reader(csvfile2)
    for row in resultsreader2:
        try:
            results4.append(map(float, row))
        except ValueError,e:
            print row
# Extract result
curves4 = np.array(results4).transpose()
print curves4


# Extract result
curves3 = np.array(results2).transpose()
print curves3
pypl.figure(1)
pypl.subplot(321)
pypl.plot(curves2[0],curves2[4]/curves2[1],label="Double Parallel Reduce")
pypl.plot(curves2[0],curves2[4]/curves2[2],label="Interval Arithmetic Filtering")
pypl.plot(curves2[0],curves2[4]/curves2[4],label="Sequential Implementation")
pypl.plot(curves2[0],curves2[4]/curves2[5],label="Total, Lazy Computation")
pypl.fill_between(curves2[0],curves2[4]/curves2[2],curves2[4]/curves2[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves2[0],curves2[4]/curves2[3],curves2[4]/curves2[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')
pypl.fill_between(curves2[0],curves2[4]/curves2[5],curves2[4]/curves2[2],label="Exact computation",alpha=0.2,color='r')

pypl.xlim(3*10**5,10**9)
pypl.ylim(0.5,5)
pypl.ylabel("Machine A")
pypl.xscale('log')
pypl.yticks(np.arange(0.5,8.5,1))
pypl.xlabel("Steep")
#pypl.legend(loc=2, prop={'size': 5})


pypl.subplot(322)
pypl.plot(curves[0],curves[4]/curves[1],label="Double Parallel Reduce")
pypl.plot(curves[0],curves[4]/curves[2],label="Interval Arithmetic Filtering")
pypl.plot(curves[0],curves[4]/curves[4],label="Sequential Implementation")
pypl.plot(curves[0],curves[4]/curves[5],label="Total, Lazy Computation")
pypl.fill_between(curves[0],curves[4]/curves[2],curves[4]/curves[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves[0],curves[4]/curves[3],curves[4]/curves[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')
pypl.fill_between(curves[0],curves[4]/curves[5],curves[4]/curves[2],label="Exact computation",alpha=0.2,color='r')

pypl.xlim(3*10**5,10**9)
pypl.ylim(0.5,4)
#pypl.ylabel("Relative Throughput")
pypl.xscale('log')
pypl.yticks(np.arange(0.5,4.5,0.5))
pypl.xlabel("Steep")

pypl.subplot(323)
pypl.plot(curves4[0],curves4[4]/curves4[1],label="Double Parallel Reduce")
pypl.plot(curves4[0],curves4[4]/curves4[2],label="Interval Arithmetic Filtering")
pypl.plot(curves4[0],curves4[4]/curves4[4],label="Sequential Implementation")
pypl.fill_between(curves4[0],curves4[4]/curves4[2],curves4[4]/curves4[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves4[0],curves4[4]/curves4[3],curves4[4]/curves4[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')

pypl.xlim(3*10**5,10**9)
pypl.ylim(0,15)
pypl.xlabel("Mss, Guarantee on Pos")
pypl.ylabel("Machine B")
pypl.xscale('log')
pypl.yticks(np.arange(0,16,2))

pypl.subplot(324)
pypl.plot(curves3[0],curves3[4]/curves3[1],label="Double Parallel Reduce")
pypl.plot(curves3[0],curves3[4]/curves3[2],label="Interval Arithmetic Filtering")
pypl.plot(curves3[0],curves3[4]/curves3[4],label="Sequential Implementation")
pypl.fill_between(curves3[0],curves3[4]/curves3[2],curves3[4]/curves3[3],label="Scheduling and Memorization Overhead",color='b',alpha=0.2)
pypl.fill_between(curves3[0],curves3[4]/curves3[3],curves3[4]/curves3[1],label="Interval Arithmetic Overhead",alpha=0.2,color='y')

pypl.xlim(3*10**5,10**9)
#pypl.ylim(0.5,4)
pypl.xlabel("Array Size")
#pypl.ylabel("Relative Throughput")
pypl.xscale('log')
#pypl.yticks(np.arange(0.5,4,0.5))
pypl.xlabel("Steep, Guarantee on Boolean b")

pypl.legend(loc=9, prop={'size': 10},bbox_to_anchor=(-0.1,-0.5))

pypl.suptitle("Speedup as a Function of Input Size",y=0.97);
#pypl.show()
pypl.savefig("Plots/figures/hybridsteepfinal.jpg");
