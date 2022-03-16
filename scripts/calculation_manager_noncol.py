#!/usr/local/bin/python3.7
# source python-select local3

import numpy as np
import os
import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Set runspinwave.sh to
# 	groundstate = 0
# 	    compile = 1
# 	  executing = 1
#forceexecuting= 0
# 	    lattice = 0
# 	 occupation = 0
#dispersionplot= 0 
# 	     spirit = 0


#
# Inputs
#

# For plotting only, 'num_mode_to_track' should not have changed.
only_plot = False

scalling_factor = 1.0

to_save_mag_structure = True

bottom = 0; top = 13
bottom = 0; top = 1.5
left = 0  ; right = 10
left = 2.2; right = 2.6

# RECOMMENDATION: In the Swiss inputcard, set 'npt' (num of kpoints) to num_kpoints_to_consider. Also decrease 'nptomega', 4 is ok.
# Try to set at most 4, otherwise, you have to increase the color and symbol lists
# num_kpoints_to_consider = 4
num_kpoints_to_consider = 1

# Number of modes (atoms in the unit cell)
num_modes = 16
# We track the lowest energy modes
num_mode_to_track = 16

# Magnetic fields
# Bfields = [i for i in np.arange(0, 1, 0.2)] + [i for i in range(1, 2, 1)] + [i for i in np.arange(2, 4, 0.4)] + [i for i in range(4, 7, 1)]
# Bfields = [i for i in np.arange(0, 4, 0.5)] + [i for i in range(4, 7, 1)]
# Bfields = [i for i in np.arange(0, 2.5, 0.2)] + [i for i in np.arange(2.5, 7, 0.5)]
# Bfields = [i for i in np.arange(0, 6, 0.1)] + [i for i in np.arange(6, 7, 0.5)] + [i for i in range(7, 9, 1)]
Bfields = [i for i in np.arange(0, 8, 0.2)]
Bfields = [i for i in np.arange(0, 1, 0.2)] + [i for i in np.arange(1, 2, 0.1)] + [i for i in np.arange(2, 3, 0.4)] + [i for i in np.arange(3, 5, 0.2)] + [i for i in np.arange(5, 8, 0.4)]
Bfields = [i for i in np.arange(0, 1, 0.2)] + [i for i in np.arange(1, 5, 0.4)] + [i for i in np.arange(5, 7, 0.2)] + [i for i in np.arange(7, 8, 0.4)]
Bfields = [i for i in np.arange(0, 0.5, 0.1)] + [i for i in np.arange(0.5, 1.5, 0.2)] + [i for i in np.arange(1.5, 2.5, 0.1)] + [i for i in np.arange(2.5, 5.5, 0.4)] + [i for i in np.arange(5.5, 7.0, 0.1)] + [i for i in np.arange(7.0, 8.0, 0.4)]

Bfields = [i for i in np.arange(0, 0.7, 0.1)] + [i for i in np.arange(0.7, 1.9, 0.2)] + [i for i in np.arange(1.9, 2.7, 0.1)] + [i for i in np.arange(2.7, 6.6, 0.2)] + [i for i in np.arange(6.6, 7.3, 0.1)] + [i for i in np.arange(7.3, 9.0, 0.2)]
Bfields = [i for i in np.arange(0, 10, 0.2)]

Bfields = [i for i in np.arange(2.6, 2.8, 0.005)]
Bfields = [2.3]

Bfield_to_save = []

if to_save_mag_structure :
   # only_plot = False

   # Bfields = [0.0, 0.4, 0.5, 1.0, 1.6, 2.0, 3.0, 4.0, 5.6, 6.0] # Set 'lattice = 1' in 'runspinwave.sh' and uncomment line 106 on 'lattice_sky.py' file
   # Bfields = [0.0, 0.3, 0.4, 0.5, 1.0, 1.6, 2.0, 2.2, 2.4, 3.0, 4.0, 5.6, 6.0, 6.9, 7.0] # Set 'lattice = 1' in 'runspinwave.sh' and uncomment line 106 on 'lattice_sky.py' file
   # Bfields = [0.0] # Set 'lattice = 1' in 'runspinwave.sh' and uncomment line 106 on 'lattice_sky.py' file
   # Bfields = [0.0, 0.4, 0.5, 1.0, 1.6, 2.0] # Set 'lattice = 1' in 'runspinwave.sh' and uncomment line 106 on 'lattice_sky.py' file
   Bfield_to_save = Bfields; lattice_filename = 'spinconfig_MnSinoncol_final.png'


# Bfields = [i for i in np.arange(0,11,0.5)] 
# Bfields = [i for i in np.arange(0,8,0.2)] 
# Bfields = [i for i in np.arange(0,1,0.1)] + [i for i in range(1,9,1)] 
# Bfields = [float(i)/4 for i in range(0,24)]
# Bfields = [0]

# IMPORTANT: various system parameter, such as anisotropy, should be sincronized between Swiss and Spirit inputcards manually
# IMPORTANT: Please make sure you don't have the keyword 'external_field_magnitude' repeated in this inputfile
spirit_inputfile           = 'input_MnSinoncol.cfg'
# IMPORTANT: Please make sure you don't have the keyword 'hm0' repeated in this inputfile
spinwave_inputfile         = 'inputcard_MnSinoncol.inp'

final_state_file           = 'spirit_output/spinconfig_MnSinoncol_final.ovf'
initial_state_file         = 'spirit_output/spinconfig_MnSinoncol_initial.ovf'
initial_state_file         = 'spirit_output/spinconfig_MnSinoncol_initial_random.ovf'

# initial_state_file         = 'spirit_output/spinconfig_MnSitest_initial.ovf'
# initial_state_file         = 'spirit_output/spinconfig_MnSitest_initialDMI.ovf'

pair_interaction_file      = 'inputfiles/pair_MnSi_rad4.5noncol.txt'
# pair_interaction_file      = 'inputfiles/pair_MnSi_rad4.5.txt'
# pair_interaction_file      = 'inputfiles/pair_MnSi_rad4.5DMI.txt'
# pair_interaction_file      = 'inputfiles/pair_MnSi_rad17.txt'
## pair_interaction_file      = 'inputfiles/pairMnSi_AFMz.txt'

output_fig_name            = 'energis_Bfield.png'
output_data_name           = 'energis_Bfield.dat'
output_intensity_name      = 'intensity_Bfield.dat'

# output_data_name           = 'energis_BfieldD0Ha.dat'
# output_data_name           = 'energis_BfieldD1Ha.dat'
# output_data_name           = 'energis_BfieldD3Ha.dat'

# output_fig_name            = 'energis_Bfield3_2.png'
# output_data_name           = 'energis_Bfield3_2.dat'

# output_fig_name            = 'energis_Bfield_saveKa_2.png'
# output_data_name           = 'energis_Bfield_saveKa_2.dat'

# output_fig_name            = 'energis_BfieldDallHa.png'
# output_data_name           = 'energis_BfieldDallHa.dat'

# output_fig_name            = 'energis_BfieldDallHc.png'
# output_data_name           = 'energis_BfieldDallHc.dat'

#
# Derivative inputs 
#

if output_data_name in [ 'energis_BfieldD0Ha.dat', 'energis_BfieldD1Ha.dat', 'energis_BfieldD3Ha.dat', 'energis_BfieldD0Hc.dat',  'energis_BfieldD1Hc.dat',  'energis_BfieldD3Hc.dat'] :
   num_kpoints_to_consider = 1

if output_data_name in [ 'energis_BfieldDallHa.dat', 'energis_BfieldDallHc.dat'] :
   only_plot = True

# Trying to get only the prefix of the file name
prefix = spinwave_inputfile[10:-4]

current_dir                = os.path.abspath(os.getcwd())

angular_momentum_file      = 'outputfiles/angular_momentum_' + prefix + '.dat'
ine_intensities_file       = 'outputfiles/ine_intensities_'  + prefix + '.dat'

pair_interaction_file_temp = pair_interaction_file[0:-4]+'_temp.txt'

# Command to run the spirit code
bashCommand_for_spirit     = 'python runspirit.py ' + spirit_inputfile + ' ' + initial_state_file + ' ' + final_state_file
# RECOMMENDATION: turn dispersion plotting in 'runspinwave.sh' off: dispersionplot=0
bashCommand_for_spinwave   = './runspinwave.sh ' + prefix 

#
# Main code
#

def editing_file(filename, key, value):
   # openning and reading the file
   file = open(filename, 'r')
   lines=file.readlines()
   file.close()

   # reopening the file in the writing mode
   target = open(filename, 'w')

   # looping through the files looking for keywords
   done = False
   for line in lines:
      if key in line and not done:
         newline = key + str(value) + '\n'
         target.write( newline )
         print( 'Modification in "' + filename + '": ' + newline  )
         done = True
      else :
         target.write( line )

   target.close()

def float_rescale(x) :
   return float(x) * scalling_factor

def rescaling_interactions(filename, filename_temp, value):
   # openning and reading the file
   file = open(filename, 'r')
   lines=file.readlines()
   file.close()

   # reopening the file in the writing mode
   target = open(filename_temp, 'w')

   # looping through the files looking for keywords
   for line_index, line in enumerate(lines) :
      
      if line_index == 0 :
         target.write( line )
      
      else :
         data_str=line.split()
         
         newline = '%3i'*5%tuple( map(int, data_str[0:5]) ) + '%16.8f'*3%tuple( map(float_rescale, data_str[5:8]) ) + '%16.8f'%( value * float(data_str[8]) ) + '\n'
         target.write( newline )

   target.close()

if not only_plot :

   # rescaling interactions if needed
   rescaling_interactions(pair_interaction_file, pair_interaction_file_temp, scalling_factor )

   # Editing the inputcards with the pair interaction file
   editing_file( spirit_inputfile  , 'interaction_pairs_file', '     ' + current_dir + '/' + pair_interaction_file_temp  )
   editing_file( spinwave_inputfile, 'pairfile'              , ' = "' + pair_interaction_file_temp + '",' )
   editing_file( spinwave_inputfile, 'basisname'             , ' = "' + final_state_file + '",'  )

   # Openning file to collect the resulting data
   output_data_file = open(output_data_name, 'w')
   output_data_file.write( '#       Bfield            second_mode_energy   first_mode_energy\n')
   output_intensity_file = open(output_intensity_name, 'w')
   output_intensity_file.write( '#       Bfield            second_mode_energy   first_mode_energy\n')

   energies_per_kpoint = []
   intensities_per_kpoint = []

   for B_index, Bfield in enumerate(Bfields) :

      # Editing spirit input file
      editing_file( spirit_inputfile, 'external_field_magnitude', '   ' + str(Bfield) )

      # Editing spinwave input file
      editing_file( spinwave_inputfile, 'hm0'                   , ' = ' + str(Bfield) + 'd0,' )

      # Running spirit to obtain relaxed spin configuration
      # print( '=================  Spirit Code  ===================')
      subprocess.call( bashCommand_for_spirit.split() )

      # Running the spinwave code
      # print( '***************** Spinwave Code *******************')
      subprocess.call( bashCommand_for_spinwave.split() )

      # Rename lattice file
      if Bfield in Bfield_to_save :
         bashCommand_cp_latticepng = 'cp ' + lattice_filename + ' ' + lattice_filename[0:-4] + '_h=' + '%.1f' % (Bfield) + '.png'
         subprocess.call(bashCommand_cp_latticepng.split())

      file = open(angular_momentum_file, 'r')
      lines=file.readlines()

      file = open( ine_intensities_file, 'r')
      lines_ine_intensities=file.readlines()

      energies = []
      intensities = []

      # Colleting the energy of modes 7 and 8 for the first few (num_kpoints_to_consider) kpoints
      for line in lines :
         data_str=line.split()

         if 'kpointindex:' in line :
            kpointindex = int(data_str[1])

         if kpointindex == num_kpoints_to_consider+1 : 
            break

         if 'mode=' in line :

            mode_index = int(data_str[1])

            if mode_index > num_modes - num_mode_to_track :
               energies.append( float(data_str[17]) + float(data_str[18]) )

               if float(data_str[18]) > 0.000 :
                  sign = -1.
               else :
                  sign =  1.

               intensity_data = lines_ine_intensities[ (kpointindex-1)*num_modes + (mode_index-1) ].split()
               intensities.append( sign*( float(intensity_data[2])+float(intensity_data[3])+float(intensity_data[4])+float(intensity_data[5]) ) )

            # print('kpoint', kpointindex, 'mode_index', mode_index, 'intensity', intensity)

      energies_per_kpoint.append(energies)
      intensities_per_kpoint.append(intensities)

      output_data_file.write( '%16.10f'%Bfield + '%16.10f'*(num_mode_to_track*num_kpoints_to_consider) % tuple( energies_per_kpoint[-1] ) + '\n')
      output_intensity_file.write( '%16.10f'%Bfield + '%16.10f'*(num_mode_to_track*num_kpoints_to_consider) % tuple( intensities_per_kpoint[-1] ) + '\n')

   output_data_file.close()

#
# Plotting
#

# Configurating the figure
params = {'font.size': 20,
          'font.family': 'serif',
          'legend.fontsize': 16,
          'mathtext.fontset': 'dejavuserif'}
plt.rcParams.update(params)
plt.clf()
plt.figure(figsize=(7.0, 5.0))
plt.subplots_adjust(left=0.14, right=0.97, top=0.95, bottom=0.2)

ax = plt.gca()
ax.tick_params(width=2,size=7)
ax.tick_params(which='both', width=1)
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=3, color='black')

ax.minorticks_on()
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)
minor_locator = AutoMinorLocator(2)
ax.yaxis.set_minor_locator(minor_locator)

# If only to plot, we need to load the data from a file
if only_plot :
   output_data_file = open(output_data_name, 'r')
   Bfields              = np.loadtxt(output_data_file)[:, 0]

   output_data_file = open(output_data_name, 'r')
   energies_per_kpoint = np.loadtxt(output_data_file)[:, 1:]

   output_intensity_file = open(output_intensity_name, 'r')
   intensities_per_kpoint = np.loadtxt(output_intensity_file)[:, 1:]

# Determining the number of kpoints, in case the read from file diverges from the variable in this script
num_kpoints_to_consider = int( len(energies_per_kpoint[0])/num_mode_to_track )
print("num_kpoints_to_consider = " , num_kpoints_to_consider)

transpose_energies_per_kpoint = list(zip(*energies_per_kpoint))
transpose_intensities_per_kpoint = list(zip(*intensities_per_kpoint))

colors  = ['black', 'red', 'forestgreen', 'dodgerblue']
colors  = ['black', 'red', 'blue', 'dodgerblue']
symbols = ['o', 's', '^', 'v']
labels = ['Q=(1, 2, 0)', 'Q=(1.009, 2, 0)', 'Q=(1.018, 2, 0)', 'Q=(1.026, 2, 0)']

if num_kpoints_to_consider == 3 :
   colors  = ['black', 'red', 'forestgreen']
   colors  = ['black', 'red', 'blue']
   symbols = ['o', 's', '^']
   labels = ['Q=(1, 2, 0)', 'Q=(1.018, 2, 0)', 'Q=(1.025, 2, 0)']

if output_fig_name == 'energis_BfieldDallHc.png' or output_fig_name == 'energis_BfieldDallHa.png' :
   labels = ['D=0.0', 'D=0.1', 'D=0.3']

cm = plt.cm.get_cmap('gnuplot')
cm = plt.cm.get_cmap('bwr')
# cm = cm.reversed()

vmax = max(max(transpose_intensities_per_kpoint))

for kpointindex in range(num_kpoints_to_consider) :

   for mode_index in range(num_mode_to_track) :

      if mode_index == 0 :
         label = labels[kpointindex]
      else :
         label = ''

      energies  = transpose_energies_per_kpoint[ kpointindex*num_mode_to_track + mode_index ]
      intensities  = transpose_intensities_per_kpoint[ kpointindex*num_mode_to_track + mode_index ]
      
      if output_data_name=="energis_BfieldDallHc.dat" and kpointindex == 2 and mode_index == 0:
         energies = np.array(energies[0:len(Bfields)-1])
         Bfields_aux = Bfields[0:len(Bfields)-1]
      else :
         Bfields_aux = Bfields

      line = plt.plot( Bfields_aux, energies, linestyle='', color  = colors[kpointindex],
                                        marker = symbols[kpointindex],
                                        label  = label,
                                        linewidth  = 0.7, 
                                        markersize = 6,
                                        fillstyle  = 'none' )[0]

      sc = plt.scatter( Bfields_aux, energies , c=intensities, vmin=-vmax, vmax=vmax, cmap=cm ) #, vmin=0, vmax=20, s=35)
      
      # line.set_clip_on(False)
      # sc.set_clip_on(False)

plt.colorbar(sc)

# plt.ylim(top=1) 
# plt.ylim(top=0.5) 
plt.ylim(top=top) 
plt.ylim(bottom=bottom)       
plt.xlim(left=left) 
plt.xlim(right=right) 


plt.xlabel('Magnetic field (T)')
plt.ylabel('Energy (meV)')

plt.legend()
plt.savefig(output_fig_name)
print ( "Saved in file: " + output_fig_name )


