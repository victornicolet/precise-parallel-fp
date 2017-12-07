import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_datafile(file_name, dtypes):
    data = np.loadtxt(file_name, delimiter=',', skiprows=0, dtype=dtypes)
    return data

def m_test_mts_errlog():
    merr_dnames = ["size2", "initmode", "inexact",
                   "superacc", "fpe2","fpe4",
                   "fpe4ee","fpe6ee","fpe8ee"]

    merr_dtypes_rec = ["int","int", "float64",
                       "float64","float64",'float64'
                       ,"float64","float64","float64"]

    merr_dtypes = { 'names' : merr_dnames,
                    'formats' : merr_dtypes_rec }

    data = read_datafile("m_test_mts_errlog.csv", merr_dtypes)

    expansions = ["inexact", "superacc","fpe2","fpe4","fpe4ee","fpe6ee","fpe8ee"]

    modes=["naive", "fpuniform", "ill conditioned"]

    plt.close('all')
    f, subplts = plt.subplots(len(modes))

    colors = ['black', 'red', '#5990aa', '#69aabb', '#aa8834', '#ff8834', '#aa6630']
    for i, mode in enumerate(modes):

        subplts[i].set_title("Data: %s" % mode)
        initmode = pd.DataFrame(data[data["initmode"] == i])
        initmode.drop("initmode",1, inplace=True)
        initmode = initmode.reset_index().groupby("size2", sort=True).mean().reset_index()

        for j, expansion in enumerate(expansions):
            subplts[i].plot(initmode['size2'], initmode[expansion],
                            linestyle = '-',
                            color = colors[j],
                            label = expansion)




    f.subplots_adjust()
    plt.show()



def m_test_experiment(fname, tname):
    m_test_mts_dnames = ["size2","initmode","inexact", "superacc","fpe2","fpe4","fpe4ee","fpe6ee","fpe8ee"]
    m_test_mts_dtypes = ["int","int","float64", "float64","float64",'float64',"float64","float64","float64"]
    mmts_dtypes = { 'names' : m_test_mts_dnames,
                    'formats' : m_test_mts_dtypes }
    data = read_datafile(fname, mmts_dtypes)



    expansions = ["inexact", "fpe2","fpe4","fpe4ee","fpe6ee","fpe8ee"]

    modes=["naive", "fpuniform", "ill conditioned"]

    plt.close('all')
    f, subplts = plt.subplots(len(modes))

    colors = ['black', 'red', '#5990aa', '#69aabb', '#aa8834', '#ff8834', '#aa6630']
    for i, mode in enumerate(modes):
        subplts[i].set_title("Data: %s" % mode)
        initmode = pd.DataFrame(data[data["initmode"] == i])
        initmode.drop("initmode",1, inplace=True)
        initmode = initmode.reset_index().groupby("size2", sort=True).mean().reset_index()

        for j, expansion in enumerate(expansions):

            subplts[i].plot(initmode['size2'], initmode[expansion],
                            linestyle = '-',
                            color = colors[j],
                            label = expansion)



    subplts[len(modes) - 1].legend(bbox_to_anchor=(0.2,0), loc="lower left",
                      bbox_transform=f.transFigure,
                                   ncol=4)
    f.subplots_adjust()
    plt.suptitle(tname, fontsize=16)
    plt.show()


m_test_experiment("m_test_poly.csv", "Polynomial evaluation - speedup/sequential.")
m_test_experiment("m_test_mts.csv", "Maxium tail sum - speedup/sequential.")
# m_test_mts_errlog()

# x = ???
# y = ???

# fig = plt.figure()

# ax1 = fig.add_subplot(111)

# ax1.set_title("Mains power stability")
# ax1.set_xlabel('time')
# ax1.set_ylabel('Mains voltage')

# ax1.plot(x,y, c='r', label='the data')

# leg = ax1.legend()

# plt.show()
