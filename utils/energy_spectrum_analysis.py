#!/usr/bin/env python

import analysis
import pandas
import matplotlib.pyplot as plt

#(iwfs, gnwfs, gewfs) = analysis.all_wfs_mapped()

#resonant:
#top = {"121nm":"121","200nm":"200","242nm":"242","400nm":"400","500nm":"500","600nm":"600","700nm":"700","800nm":"800","900nm":"900","1000nm":"1000","1100nm":"1100","1200nm":"1200"}
#top = { 1.5E14:"1.5E14",  1E13:"1E13",   2E13:"2E13",   2E14:"2E14",   3E13:"3E13",   4E13:"4E13",   5E13:"5E13",   6E13:"6E13",   7E13:"7E13",   8E13:"8E13",   9E13:"9E13"}
top = { 800:"resonant_long", 121:"resonant_exact_long", 1200:"non-resonant_long" }
folders = ["hydrogen","hydrogen_fs","hydrogen_sfa"]
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
        iwf = analysis.indexed_wf( analysis.get_prototype(top[j] + "/" + i + "/prototype.csv"), analysis.import_petsc_vec( top[j] + "/" + i + "/wf_final.dat") )
        ion = analysis.ionization(iwf)
        print ( top[j] + " - " + i + " -- " + str(1.-ion["bound"]) + ", absorbed: " + str(ion["absorbed"]) )
        ionization[j][i] = 1.-ion["bound"]
        gnwf = analysis.grouped_by("n", iwf)
        gewf = analysis.binned_group_by(iwf, .01)
        de[i] = gewf["wf_abs"]
        dn[i] = gnwf["wf_abs"].loc[2:28]
        indexn = gnwf.loc[:28].index
        indexe = gewf.index

    datae[j] = pandas.DataFrame(de, indexe)
    datan[j] = pandas.DataFrame(dn, indexn)

    plt.figure()
    datan[j].plot(logy=False, title=top[j])
    plt.savefig(top[j] + "_n.pdf")
    
    plt.figure()
    datae[j].plot(logy=True, title=top[j])
    plt.savefig(top[j] + "_e.pdf")

ion = pandas.DataFrame( ionization ).transpose().sort()
plt.figure()
ion.plot(logy=True, logx=True, title="ionization vs. intensity")
plt.savefig("ionization.pdf")

for i in folders:
    de_l = {}
    dn_l = {}
    for j in top:
        de_l[j] = datae[j][i]
        dn_l[j] = datan[j][i]
    de = pandas.DataFrame( de_l, datae[1.5E14].index)
    dn = pandas.DataFrame( dn_l, datan[1.5E14].index)
    plt.figure()
    de.plot(logy=True, title=i)
    plt.savefig(i + "_e.pdf")
    plt.figure()
    dn.plot(logy=True, title=i)
    plt.savefig(i + "_n.pdf")

