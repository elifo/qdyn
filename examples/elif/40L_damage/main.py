### 10Lnuc lvfz model of Idini & Ampuero paper        ###
### modified for Python wrapper after B.Idini's files ###
### for comparison with Florez et al., 
### damage half-width set to 3 m (2h/l_vw=40).
###
### Elif (06/2023) github@elifo ###
###
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pickle
import os
import sys
from numpy import random
import time

# Add QDYN directory to PATH
sys.path.append('/Users/elifo/Desktop/2023_Cycles/qdyn-read-only/src')  # For pyqdyn
sys.path.append('/Users/elifo/Desktop/2023_Cycles/qdyn-read-only/utils/post_processing') # For plot_functions
print(sys.path)
from pyqdyn import qdyn
import plot_functions as qdyn_plot

start = time.time()

# Instantiate the QDYN class object
p = qdyn()

# Predefine parameters
t_yr = 3600.0* 24.0* 365.0    # seconds per year
Lasp = 40                     # Length of asperity (VW) / nucleation length
L = 4                         # Length of model / Lasp
a = 0.014
b = 0.019
resolution = 5                # Mesh resolution / process zone width

# Get the settings dict
set_dict = p.set_dict

""" Step 1: Define simulation/mesh parameters """
# Global simulation parameters
set_dict["MESHDIM"] = 1        # Simulation dimensionality (1D fault in 2D medium)
set_dict["FINITE"] = 1         # 
set_dict["TMAX"] = 10*t_yr      # Maximum simulation time [s]
set_dict["V_PL"] = 1e-9        # Plate velocity
set_dict["MU"] = 3e10          # Shear modulus
set_dict["SIGMA"] = 120e6        # Effective normal stress [Pa]
set_dict["ACC"] = 1e-10         # Solver accuracy
# set_dict["SOLVER"] = 1         # Solver type (Runge-Kutta)

# Setting some (default) RSF parameter values
set_dict["SET_DICT_RSF"]["A"] = a          # Direct effect (will be overwritten later)
set_dict["SET_DICT_RSF"]["B"] = b          # Evolution effect
set_dict["SET_DICT_RSF"]["DC"] = 2e-3      # Characteristic slip distance
set_dict["SET_DICT_RSF"]["V_SS"] = 1e-9    # Reference velocity [m/s]
set_dict["SET_DICT_RSF"]["TH_0"] = set_dict["SET_DICT_RSF"]["DC"] / set_dict["V_PL"]    # Initial state [s]

# Damage zone
set_dict['D'] = 0.9 # damage level (delta) = 1- (mu_damage)/ mu_intact_rock 
set_dict['HD'] = 12.0

# Compute relevant length scales:
mu_damage = (1.0- set_dict['D'])* set_dict["MU"]
Lc = 2.0/ np.pi* mu_damage* set_dict["SET_DICT_RSF"]["DC"]* b/  \
                     (set_dict["SIGMA"]* (b-a)** 2.0)
print ('Nucleation length (m): ', Lc)

Lb = 9.0* np.pi/ 32.0* mu_damage* set_dict["SET_DICT_RSF"]["DC"]/ \
                      (set_dict["SIGMA"]* b)
print ('Process zone width (m): ', Lb)

# Length of asperity [m]
Lasp *= Lc
# Fault length [m]
L *= Lasp

print ('VW patch (asperity) size, Lasp (m, Lc): ', Lasp, Lasp/Lc)
print ('Model size, L (m, Lc): ', L, L/Lc)
print ('*')
print ('Damage width ratio, 2h/L_vw: ', 2.0* set_dict['HD']/ Lasp)
print ('Damage level, delta (%): ', 100.0* set_dict['D'])

# Find next power of two for number of mesh elements
N = int(np.power(2, np.ceil(np.log2(resolution * L / Lb))))
# Spatial coordinate for mesh
x = np.linspace(-L/2, L/2, N, dtype=float)

# Set mesh size and fault length
set_dict["N"] = N
set_dict["L"] = L
# Set time series output node to the middle of the fault
set_dict["IC"] = N // 2

# output settings
set_dict['NXOUT'] = max(1, set_dict['NXOUT']/1024) # Snapshot resolution (every N elements)
set_dict["NTOUT"] = 10                             # Save output every N steps


""" Step 2: Set (default) parameter values and generate mesh """
p.settings(set_dict)
p.render_mesh()

""" Step 3: override default mesh values """
# Distribute direct effect a over mesh according to some arbitrary function
# outside VW patch, set a = (3b-a) after benjamin's
cdt = (abs(x) > Lasp/ 2.0)
p.mesh_dict["A"][cdt] = 3.0* p.mesh_dict["B"][cdt]- a
print (p.mesh_dict["A"][cdt])

# Write input to qdyn.in
p.write_input()



### RUN ###
p.run()
end = time.time()
print ('*** Elapsed time (min): ', (end- start)/60.0 ); print ()


### POST-PROCESS ###
p.read_output()

# Time series of stress, state, and maximum slip rate on the fault
qdyn_plot.timeseries(p.ot[0], p.ot_vmax)

# Spatio-temporal evolution of slip rates
qdyn_plot.slip_profile(p.ox, warm_up=1*t_yr)

### FIN.