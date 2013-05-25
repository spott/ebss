#!/usr/bin/env python

import analysis
import pandas
import matplotlib.pyplot as plt

top = { 3:"3", 6:"6", 9:"9", 12:"12", 15:"15", 18:"18", 21:"21"}
folders = ["hydrogen"]
datae = {}
datan = {}
ionization = {}
print ("ionization:")
for j in top:
    de = {}
    dn = {}
    indexn = []
    indexe = []
    ionization[j] = {}
    for i in folders:
        analysis.gen_k_spectrum_figure( top[j] + "/wf_final_kspectrum.dat" , .5, top[j] + "_" + i + "kspectrum.pdf" )
        iwf = analysis.indexed_wf( analysis.get_prototype(top[j] + "/prototype.csv"), analysis.import_petsc_vec( top[j] + "/wf_final.dat") )
        ion = analysis.ionization(iwf)
        print ( top[j] + " - " + i + " -- " + str(1.-ion["bound"]) + ", absorbed: " + str(ion["absorbed"]) )
        ionization[j][i] = 1.-ion["bound"]
        gnwf = analysis.grouped_by("n", iwf)
        gewf = analysis.binned_group_by(iwf, .01)
        de[i] = gewf["wf_abs"]
        dn[i] = gnwf["wf_abs"].loc[:28]
        indexn = gnwf.loc[:28].index
        indexe = gewf.index

    datae[j] = pandas.DataFrame(de, indexe)
    datan[j] = pandas.DataFrame(dn, indexn)

    plt.figure()
    datan[j].plot(logy=True, title=top[j])
    plt.savefig(top[j] + "_n.pdf")
    
    plt.figure()
    datae[j].plot(logy=True, title=top[j])
    plt.savefig(top[j] + "_e.pdf")

ion = pandas.DataFrame( ionization ).transpose().sort()
plt.figure()
ion.plot(logy=True, logx=True, title="ionization vs. cycles")
plt.savefig("ionization.pdf")

for i in folders:
    de_l = {}
    dn_l = {}
    for j in top:
        de_l[j] = datae[j][i]
        dn_l[j] = datan[j][i]
    de = pandas.DataFrame( de_l, datae[3].index)
    dn = pandas.DataFrame( dn_l, datan[3].index)
    plt.figure()
    de.plot(logy=True, title=i)
    plt.savefig(i + "_e.pdf")
    plt.figure()
    dn.plot(logy=True, title=i)
    plt.savefig(i + "_n.pdf")

