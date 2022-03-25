#!/usr/local/bin/python3 
import os
import sys
import socket

hostname = socket.gethostname()
print('runspirit.py - Hostname:', hostname)

### Make sure to find the Spirit modules
### This is only needed if you did not install the package
# spirit_py_dir = os.path.dirname(os.path.realpath(__file__)) + "core/python/Spirit"
# spirit_py_dir = os.path.abspath(os.path.join("~/applications/spirit-develop/core/python"))
# spirit_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "../../../../../Applications/spiritGit/core/python"))
# spirit_py_dir = '/Users/santos/Applications/spiritGit/core/python'  # mb-santos
# spirit_py_dir = '/Users/flavianojs/Applications/spirit/core/python' # Flaviano's MacBook

# if   hostname == 'theospc47' :
#     spirit_py_dir = '/home/fdossantos/codes/spirit/core/python'

# elif hostname == 'Flavianos-MacBook-Pro.local'  or hostname == 'tsf-452-wpa-4-009.epfl.ch' :
#     spirit_py_dir = '/Users/flavianojs/Applications/spirit/core/python'

# elif hostname == 'mb-dossantos' :
#     spirit_py_dir = '/Users/santos/Applications/spiritGit/core/python'

# else :
#     print('Unknown Hostcomputer:', hostname)
#     exit()

print ('runspirit.py - Spirit installed via pip' )
# print ('runspirit.py - Spirit code directory: ' + spirit_py_dir)
# sys.path.insert(0, spirit_py_dir)

### Import numpy
import numpy as np

### Import Spirit modules
from spirit import state
from spirit import hamiltonian
from spirit import system
from spirit import geometry
from spirit import chain
from spirit import configuration
from spirit import transition
from spirit import simulation
from spirit import quantities
from spirit import io
from spirit import log
from spirit import parameters


GNEB = simulation.METHOD_GNEB
LLG  = simulation.METHOD_LLG
Heun = simulation.SOLVER_HEUN
VP   = simulation.SOLVER_VP

convergence = 10e-14

# Reading the inputfile
cfgfile = sys.argv[1]
p_state=state.setup(cfgfile)

# Loading the initial state
initial_state_file = sys.argv[2]
io.image_read(p_state, initial_state_file)

# parameters.llg.set_convergence(p_state, convergence, idx_image=-1, idx_chain=-1)

simulation.start(p_state, method_type=LLG, solver_type=VP, n_iterations=10000000, n_iterations_log=-1, single_shot=False, idx_image=-1, idx_chain=-1)

spin_config_filename = sys.argv[3]
io.image_write(p_state, spin_config_filename)
# print("OUTPUT: " + spin_config_filename )

