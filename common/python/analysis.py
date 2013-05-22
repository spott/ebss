import numpy
import pandas
import struct
import matplotlib.pyplot as plt
import os

#we are interested in the energy levels:

#prototype load:
def get_prototype( filename = "prototype.csv" ):
    prototype_f = file( filename )

    prototype = []

    for l in prototype_f:
        i = l.split(',')
        a = dict()
        a["n"] = int(i[0])
        a["j"] = float( int( i[1] ) ) /2.
        a["l"] = int(i[2])
        a["e"] = float(i[4])
        prototype.append(a)
    return prototype

def import_petsc_vec( filename ):
    with open(filename, "rb") as f:
        byte = f.read(8)
        size = struct.unpack('>ii',byte)[1]
        npy = numpy.fromfile(f, '>D', size)
        return npy

def indexed_wf( prototype, wf ):
    r = pandas.DataFrame( prototype )
    r['wf'] = pandas.Series( wf )
    return r

def grouped_by( index, indexedwf ):
    #get an absolute value for the wf, otherwise we can't aggregate:
    r = indexedwf.copy()
    sqabs = lambda x : abs(x)**2
    r["wf"] = r["wf"].map( sqabs )
    g = r.groupby(index)
    r = g.aggregate( pandas.np.sum )
    return pandas.DataFrame({"wf_abs": r['wf']})

def binfunc( e, bins ):
    for b in bins:
        if e < b[0]:
            return b[1]
        else:
            continue
    return bins[-1][1]

def binned_group_by( bin_size, indexedwf ):
    r = indexedwf.copy()
    r["wf"] = r["wf"].map( lambda x: abs(x)**2 )
    minimum = min(r["e"])
    maximum = max(r["e"])
    diff = maximum - minimum
    #create the bins:
    bins = [ ((i + 1) * bin_size + minimum, (i + .5) * bin_size + minimum ) for i in range(int ( diff / bin_size ) + 1) ]
    bf = lambda x: binfunc( x, bins )
    g = r.set_index('e').groupby( bf )
    r = g.aggregate( pandas.np.sum )
    return pandas.DataFrame({"wf_abs":r['wf']})

def wavefunctions():
    filenames = os.listdir("./")
    wf_filenames = 0
    for i in filenames:
        if (i[:2] == "wf"):
            wf_filenames += 1
    wfs = []
    for i in range(wf_filenames-1):
        wfs.append( import_petsc_vec( "wf_" + str(i) + ".dat" ) )
    wfs.append( import_petsc_vec( "wf_final.dat" ) )
    return wfs

def ionization( indexedwf, cutoff = 0.0 ):
    r = indexedwf.copy()
    r["wf"] = r["wf"].map( lambda x: abs(x)**2 )
    bf = lambda x:  True if x > cutoff else False
    g = r.set_index('e').groupby(bf)
    r = g.aggregate( pandas.np.sum )
    return {"bound" : r["wf"].loc[False], "ionized" : r["wf"].loc[True], "absorbed": 1- r["wf"].sum()}


def all_wfs_mapped():
    prototype = get_prototype()
    wfs = wavefunctions()
    iwfs = [ indexed_wf(prototype, wf) for wf in wfs]
    gnwfs= [ grouped_by("n", wf) for wf in iwfs]
    gewfs = [binned_group_by(.01, wf) for wf in iwfs]
    return (iwfs, gnwfs, gewfs)

def final_wf_mapped( bin_size = .01):
    prototype = get_prototype()
    iwf = indexed_wf( prototype, import_petsc_vec( "wf_final.dat" ) )
    gnwf= grouped_by("n", iwf)
    gewf = binned_group_by(bin_size, iwf)
    return (iwf, gnwf, gewf)

def timeseries(fname, tfname = "time.dat", title="tsdata" ):
    with open(tfname, "rb") as f:
        time = numpy.fromfile(f, 'd', -1)
    with open(fname, "rb") as f:
        ef = numpy.fromfile(f, 'd', -1)
    return pandas.DataFrame( {title: ef}, index=time)

def efield( effname = "efield.dat" ):
    return timeseries( fname = effname ,title="efield")


