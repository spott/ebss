#!/usr/bin/env python

import itertools
import os
import shlex, subprocess

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

R_max_list = [ i for i in range(1000, 7000, 1000) ]
n_max_list = [ 50 * i for i in range(1, 6) ]

product = itertools.product(R_max_list, n_max_list)

for i in product:
    R_max = i[0]
    n_max = i[1]
    basis_folder = "./" + str(R_max) + "/" + str(n_max) + "/"
    hamiltonian_folder = "./" + str(R_max) + "/" + str( n_max ) + "_hamiltonian/"
    nonlinear_folder = "./" + str( R_max ) + "/" + str( n_max ) + "_nonlinear/"

    ensure_dir(basis_folder)
    ensure_dir(hamiltonian_folder)
    ensure_dir(nonlinear_folder)

    b = open( basis_folder + "run", "w" )
    b.write( "findbasis -basis_atom hydrogen -basis_nmax " + str(n_max) + " -basis_lmax 5 -basis_rmax " + str(R_max))
    os.fchmod( b.fileno(), 0o0777)
    b.close()

    h = open( hamiltonian_folder + "run", "w" )
    h.write( "mpiexec -n 12 findhamiltonian -hamiltonian_basis_config ../" + str(n_max) + "/Basis.config")
    os.fchmod( h.fileno(), 0o0777)
    h.close()

    n = open( nonlinear_folder + "run", "w" )
    n.write( "mpiexec -n 12 nonlinear_index -hamiltonian_config ../" + str(n_max) + "_hamiltonian/Hamiltonian.config")
    os.fchmod( n.fileno(), 0o0777)
    n.close()

current_dir = os.getcwd()
#print(current_dir)

#create bases:
p = itertools.product(R_max_list, n_max_list)

for j in p:
    R_max = j[0]
    n_max = j[1]
    basis_folder = os.path.abspath("./" + str(R_max) + "/" + str(n_max) + "/")
    os.chdir(basis_folder)
    #run the file:
    f = open( "run" , 'r')
    command = shlex.split(f.read())
    f.close()
    print command

    f = open( "err" , 'w')
    subprocess.call(command , stderr=f.fileno() )
    f.close()

    os.chdir(current_dir)

    hamiltonian_folder = os.path.abspath("./" + str(R_max) + "/" + str(n_max) + "_hamiltonian/")
    os.chdir(hamiltonian_folder)

    f = open( "run" , 'r')
    command = shlex.split(f.read())
    f.close()
    print command

    f = open( "err" , 'w')
    subprocess.call(command , stderr=f.fileno() )
    f.close()

    os.chdir(current_dir)

    nonlinear_folder = os.path.abspath("./" + str(R_max) + "/" + str(n_max) + "_nonlinear/")
    os.chdir(nonlinear_folder)

    f = open( "run" , 'r')
    command = shlex.split(f.read())
    f.close()
    print command

    f = open( "err" , 'w')
    subprocess.call(command , stderr=f.fileno() )
    f.close()

    os.chdir(current_dir)


