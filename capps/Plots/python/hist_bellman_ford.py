import csv
import matplotlib.pyplot as pypl
import numpy as np

pypl.style.use('ggplot')

results = []
with open("Plots/bellman_ford.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

barw = 0.25

pypl.bar(0.25,1, color='red', width=barw,label='Floating-Point')

pypl.bar([0.75],[results[4][0]], color = 'blue', width=barw,label='Interval Arithmetic')
pypl.bar([1],[results[5][0]],width=barw,color = 'grey',label='Reverse processing')
pypl.bar([1.25],[results[3][0]],color='green', width=barw,label='Exact Computation')
pypl.bar([1.5],[results[6][0]],width=barw,label='Total')

pypl.xticks([])
pypl.yticks(np.arange(0, 10, step=0.5))
pypl.ylabel("Relative Computation Time")

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Relative Mean Computation Time for Bellman-Ford");
#pypl.show()
pypl.savefig("Plots/histBellmanFord.jpg");

