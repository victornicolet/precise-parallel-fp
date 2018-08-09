import csv
import matplotlib.pyplot as pypl
import numpy as np


results = []
with open("Plots/csv/lazympsseqalt.csv", "rb") as csvfile:
    resultsreader = csv.reader(csvfile)
    for row in resultsreader:
        results.append(map(float, row))


results = np.array(results).transpose()
print results

barw = 0.25

pypl.bar(0.25,1, color='red', width=barw,label='Floating-Point')

pypl.bar([0.75,2],[results[5][0],results[5][0]], color = 'blue', width=barw,label='Filtering Step')
pypl.bar([1,2.25],[results[6][0],results[7][0]],width=barw,color = 'grey',label='Reverse Step')
#pypl.bar([1.25,2.5],[results[8][0],results[9][0]],color='green', width=barw/2,label='Exact Computation')
pypl.bar([1.25,2.5],[results[16][0],results[17][0]],color='green', width=barw,label='Exact Computation Step')
pypl.bar([1.5,2.75],[results[12][0],results[13][0]],width=barw,label='Total',color='yellow')

pypl.xticks([1.25,2.5],['Guarantee on Mps','Guarantee on Pos'])
#pypl.yticks(np.arange(0, 10, step=0.5))
pypl.ylabel("Relative Wall Time")
pypl.xlim(0,3.25)

pypl.legend(loc=2, prop={'size': 10})
pypl.title("Relative Mean Computation Time for Maximum Prefix Sum");
#pypl.show()
pypl.savefig("Plots/figures/histmpsseqalt.jpg");

