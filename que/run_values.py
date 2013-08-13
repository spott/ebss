#!/usr/bin/env python

import itertools
import os
import datetime
import shlex, subprocess

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

R_max_list = [ i for i in range(500, 1000, 500) ]
n_max_list = [ 50 * i for i in range(1, 28, 2) ]
points_list = [ 10000 * i for i in range(2, 9 ,2) ]

current_dir = os.getcwd()
proc = 12
atom = "hydrogen"
que = "crc-serial"
lmax = 5

que = """#!/bin/bash

#PBS -A UCB00000170
#PBS -q {que}
#PBS -N {atom}_{t}
#PBS -l walltime={time}
#PBS -l nodes=1:ppn={proc}
#PBS -m abe
#PBS -M andrew.spott@gmail.com
{other}

module load openmpi/openmpi-1.6.4_gcc-4.7.2_torque-4.2.3_ib
module load torque/torque-4.2.3
module load moab/moab-7.2.2
module load gsl/gsl-1.15_gcc-4.7.2
module load python/epd-7.3.2
module load fftw/fftw-3.3.3_openmpi-1.6.4_gcc-4.7.2_torque-4.2.3_double

cd $PBS_O_WORKDIR

{script}"""

for i in itertools.product(R_max_list, points_list):
    R_max = i[0]
    points = i[1]
    basis_folder = "./" + str(R_max) + "/" + str(points) + "/basis/"
    print("finding basis for folder" + basis_folder)
    if (not os.path.exists(basis_folder + "/done") and ((points > 70000 and R_max > 500) or (R_max < 1000)) ):
        ensure_dir(basis_folder)
        b = open( basis_folder + "que", "w" )
        run_command = "findbasis -basis_atom " + atom + " -basis_nmax " + str(max(n_max_list)) + " -basis_points " + str(points) + " -basis_lmax "+lmax+" -basis_rmax " + str(R_max) + "\n" + "touch done \n"
        b.write( que.format( atom=atom, proc=str(lmax), script=run_command, que=que, other="", t="basis", time="24:00:00" ))
        os.fchmod( b.fileno(), 0o0777)
        b.close()
        basis_folder = os.path.abspath(basis_folder)
        os.chdir(basis_folder)
        #run the file:
        #f = open( "que" , 'r')
        #command = shlex.split(f.read())
        f.close()
        print command

        err = open( "err" , 'w')
        out = open( "proc_id" , 'w')
        subprocess.call(['qsub', 'que']  , stderr=err.fileno(), stdout=out.fileno() )
        err.close()
        out.close()

    print("done")

    os.chdir(current_dir)

current_dir = os.getcwd()
#print(current_dir)

#create bases:
p = itertools.product(R_max_list, n_max_list, points_list)

for j in p:
    R_max = j[0]
    n_max = j[1]
    points = j[2]

    os.chdir(current_dir)
    with file("./" + str(R_max) + "/" + str(points) + "/proc_id") as f:
        basis_proc_id = f.readline()
    hamiltonian_folder = os.path.abspath("./" + str(R_max) + "/" + str(points) + "/" + str(n_max) + "_hamiltonian/")
    ensure_dir(hamiltonian_folder)
    os.chdir(hamiltonian_folder)

    h_p_id = str()
    if ( not os.path.exists("done") ):
        h = open( "que", "w" )
        hscript = "mpiexec -n " + str(proc) + "  findhamiltonian -hamiltonian_basis_config ../basis/Basis.config -hamiltonian_nmax " + str(n_max) + "\ntouch done \n"
        h.write( que.format( atom=atom, proc=str(proc), script=hscript, que=que, other="#PBS -W after:"+basis_proc_id, t="hamiltonian", time="24:00:00" ))
        os.fchmod( h.fileno(), 0o0777)
        h.close()

        err = open( "err" , 'w')
        subprocess.call(["qsub" "que" , stderr=err.fileno(), stdout=h_p_id )
        err.close()

        f.write(str(command))
        f.write(str(datetime.datetime.today()))
        f.close()

    os.chdir(current_dir)

    nonlinear_folder = os.path.abspath("./" + str(R_max) + "/" + str(points) + "/" + str(n_max) + "_nonlinear/")
    ensure_dir(nonlinear_folder)
    os.chdir(nonlinear_folder)
    print( nonlinear_folder )

    if ( not os.path.exists("done") ):
        n = open( "que", "w" )
        nscript = "mpiexec -n " + str(proc) + " nonlinear_index -hamiltonian_config ../" + str(n_max) + "_hamiltonian/Hamiltonian.config -nonlinear_chi1 1 -nonlinear_chi3 1,0,0 -nonlinear_chi3 1,1,1 -nonlinear_chi3 1,-1,1 -nonlinear_chi5 1,0,0,0,0 -nonlinear_chi5 1,1,1,-1,-1 -nonlinear_chi5 1,1,1,1,1 -nonlinear_freq 0 -nonlinear_wavelengths 400,800,267,1200,1800,4000 \ntouch done \n"
        n.write( que.format( atom=atom, proc=str(proc), script=nscript, que=que, other="#PBS -W after:"+h_p_id, t="nonlinear", time="24:00:00" ))
        os.fchmod( n.fileno(), 0o0777)
        n.close()

        err = open( "err" , 'w')
        out = open( "proc_id" , 'w')
        subprocess.call(["qsub","que"] , stderr=err.fileno(), stdout=out.fileno() )
        err.close()
        out.close()
        

    os.chdir(current_dir)


