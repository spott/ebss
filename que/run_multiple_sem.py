#!/usr/bin/env python2.7

import argparse
import os.path
import itertools
import os
import subprocess
import stat
import re
import time
import random
import string

#def grouper(n, l):
        #return [l[j:j+n] for j in range(0, len(l), n)]

def grouper(n, iterable):
    i = iter(iterable)
    piece = list(itertools.islice(i, n))
    while piece:
        yield piece
        piece = list(itertools.islice(i, n))
#def grouper(n, iterable, padvalue=None):
    #"""grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"""
    #return itertools.izip(*[itertools.chain(iterable, itertools.repeat(padvalue, n-1))]*n)

parser = argparse.ArgumentParser(description='Run a number of different jobs with a set of common configs, but small differences')

#Default config files:
parser.add_argument('-output_folder','-o', default=os.path.abspath("./"),
                      help='the base folder under which all variations will be put')
parser.add_argument('-hamiltonian_config', default="./Hamiltonian.config",
                      help='The config file with all the defaults for the Hamiltonian')
parser.add_argument('-laser_defaults', default="./Laser.config",
                      help='The config file with all the defaults for the laser')
parser.add_argument('-absorber_defaults', default="./Absorber.config",
                      help='The config file with all the defaults for the absorber')
parser.add_argument('-state_defaults', default="./State.config",
                      help='The config file with all the defaults for the state')
parser.add_argument('-pulsetrain_defaults', default="./Pulsetrain.config",
                      help='The config file with all the defaults for the pulsetrain')
parser.add_argument('-dipole_defaults', default="./Dipole.config",
                      help='The config file with all the defaults for the pulsetrain')

#PBS stuff:
parser.add_argument('-ppn', type=int, default=12,
                      help='Number of processors per node')
parser.add_argument('-nodes', type=int, default=4,
                      help='Number of nodes')
parser.add_argument('-time', default="4:00",
                      help='Walltime')
parser.add_argument('-que', default="janus-short",
                      help='que to use')
parser.add_argument('-jobname_prefix', default="",
                      help='the prefix for the job name')
#parser.add_argument('-que_args',default="", help="arguments to add to the que")

#parameters to iterate over, add to this as I need it:
parser.add_argument('--param_linear', '-l', nargs='*',action='append',
        help='The key value and the max/min/delta for that value')
parser.add_argument('--param_list', '-i', nargs='*',action='append',
        help='The key value and the values for that parameter')
#parser.add_argument('-pulsetrain_spacing

args = parser.parse_args()
argdict = vars(args)

argdict["hamiltonian_config"] = os.path.abspath(argdict["hamiltonian_config"])
argdict["laser_defaults"] = os.path.abspath(argdict["laser_defaults"])
argdict["absorber_defaults"] = os.path.abspath(argdict["absorber_defaults"])
argdict["state_defaults"] = os.path.abspath(argdict["state_defaults"])
argdict["pulsetrain_defaults"] = os.path.abspath(argdict["pulsetrain_defaults"])
argdict["dipole_defaults"] = os.path.abspath(argdict["dipole_defaults"])
argdict["jobname_prefix"] = "propagate" #argdict["jobname_prefix"]+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))

#print(argdict)
que_template = """#!/bin/bash

#PBS -A UCB00000170
#PBS -N {jobname}
#PBS -l walltime={time}
#PBS -l nodes=1:ppn={ppn}
#PBS -m abe
#PBS -M andrew.spott@gmail.com

. /curc/tools/utils/dkinit

reuse .openmpi-1.6_gcc-4.7.1_torque-2.5.11_ib
reuse GSL
reuse Parallel

cd /home/ansp6066/code/ebss/propagate/

#make clean
#make all

echo "


"

git log | head -n 10

cd $PBS_O_WORKDIR
"""

run_template = """
mkdir {directory}

sem --id $PBS_JOBID -j{proc} --files cd {directory} ";" propagate \\
    -hamiltonian_config {hamiltonian} \\
    -laser_config {laser_defaults} \\
    -absorber_config {absorber_defaults} \\
    -state_config {state_defaults} \\
    -pulsetrain_config {pulsetrain_defaults} \\
    -dipole_config {dipole_defaults} \\
{parameters}    -not_shared_tmp \\
    -viewer_binary_skip_info \\
    -log_summary ";" \\
    cd ../ ";" tar -zcf {directory}.tar.gz {directory}/ "&&"\\
    rm -rf {directory}/ ";"
"""

#make tensor product list:
parameters = list()
for p in argdict["param_list"]:
    ptype = list()
    for l in p[1:]:
        ptype.append(["-" + p[0],l])
    parameters.append(ptype)

parameter_sets = itertools.product(*parameters)

filenum = 0
#make the que text:
pss = grouper(args.nodes * args.ppn, parameter_sets)
for a in pss:
    qftext = que_template.format(
            jobname = argdict["jobname_prefix"],
            time = args.time,
            nodes = args.nodes,
            ppn = args.ppn)

    for p in a:
        directory = ""
        paramstring = ""
        for i in p:
            directory += i[0][1:] + "_" + i[1] + "_"
            paramstring += "    " + i[0] + " " + i[1] + " \\\n"
        qftext += run_template.format(
                key = argdict["jobname_prefix"],
                directory = directory,
                proc = args.ppn,
                hamiltonian = argdict["hamiltonian_config"],
                laser_defaults = argdict["laser_defaults"],
                absorber_defaults = argdict["absorber_defaults"],
                state_defaults = argdict["state_defaults"],
                pulsetrain_defaults = argdict["pulsetrain_defaults"],
                dipole_defaults = argdict["dipole_defaults"],
                parameters = paramstring)

    qftext += "\nsem --id $PBS_JOBID --wait"
    qf = file("que" + str(filenum),'w')
    qf.write(qftext)
    filenum += 1

exit(0)

#compile the code:
#try:
    #dir = os.getcwd()
    #os.chdir("/home/ansp6066/code/ebss/propagate/")
    #out = subprocess.check_output(["make","clean"])
    #out += subprocess.check_output(["make","all"])
    #print(out)
    #os.chdir(dir)
#except subprocess.CalledProcessError:
    #print("make returned an error")
    #exit()

#for p in parameter_sets:
    #print p
    #create a directory name:
    #directory = ""
    #paramstring = ""
    #for i in p:
        #directory += i[0][1:] + "_" + i[1] + "_"
        #paramstring += "    " + i[0] + " " + i[1] + " \\\n"
        ##directory = re.sub(r',',r'-', directory, count=0, flags=0)
    #print(directory)
    ##make the directory:
    ##try:
        ##os.makedirs(directory)
    ##except OSError:
        ##print("couldn't make the directory for some reason: " + directory )
    ##make the que file:
    #qf = file(directory + "/que", 'w')
    #qf.write(que_template.format(
        #jobname = args.jobname_prefix + "propagate",
        #time = args.time,
        #nodes = args.nodes,
        #ppn = args.ppn)
        ##proc = args.nodes * args.ppn,
        ##hamiltonian = argdict["hamiltonian_config"],
        ##laser_defaults = argdict["laser_defaults"],
        ##absorber_defaults = argdict["absorber_defaults"],
        ##state_defaults = argdict["state_defaults"],
        ##pulsetrain_defaults = argdict["pulsetrain_defaults"],
        ##parameters = paramstring))

    #for 
    #qf.close()
    #out = ""
    #try:
        #os.chdir(directory)
        #os.chmod("que", stat.S_IRWXU | stat.S_IRWXG)
        ##print(args.que_args)
        #out = subprocess.check_output(["qsub", "-q", args.que,"./que"])
        #procidf = file("proc_id", 'w')
        #procidf.write(out)
        #procidf.close()
        #os.chdir("../")
    #except subprocess.CalledProcessError:
        #print ("qsub returned an error")
    #time.sleep(1)





def frange(limit1, limit2 = None, increment = 1.):
    """
    Range function that accepts floats (and integers).

    Usage:
    frange(-2, 2, 0.1)
    frange(10)
    frange(10, increment = 0.5)

    The returned value is an iterator.  Use list(frange) for a list.
    """

    if limit2 is None:
        limit2, limit1 = limit1, 0.
    else:
        limit1 = float(limit1)

    count = int(math.ceil(limit2 - limit1)/increment)
    return (limit1 + n*increment for n in range(count))


