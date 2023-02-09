# Run with: python kpath.py XX
# XX is the inputcard sufix identifier

import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams

params = {
          # 'font.size': 18,
          'font.family': 'serif',
          # 'legend.fontsize': 14,
          'figure.figsize': (5.5, 5) ,
          "legend.loc": 'upper left',
          'mathtext.fontset': 'dejavuserif'}
plt.rcParams.update(params)

plt.clf()
axes = plt.gca()

fileinput = "outputfiles/kpath_"+sys.argv[1]+".dat" # "outputfiles/kpath_XX.dat"
print( 'kpath input file:', fileinput )
file = open(fileinput)
lines_input_file=file.readlines()

counter = 0

points_x = []
points_y = []
for line in lines_input_file :
   data_str=line.split()

   toplot = False
   if data_str[0].startswith('#') :
      toplot = True
      counter += 1
      # print( "PLOT", counter )
      if counter == 1 :
         plt.plot(points_x, points_y, linewidth=1.5, color='black', linestyle='', marker='o', markersize=5, label='Reciprocal lattice' )
         box_lim = max(points_x)*1.1
      elif counter == 2 :
         plt.plot(points_x, points_y, linewidth=1.5, color='red', linestyle='-', marker='', label=r'$k$ path' )
      elif counter == 3 :
         plt.arrow(points_x[0], points_y[0], points_x[1], points_y[1], color='blue', linestyle='--', head_width=0.08, length_includes_head=True)

         text_pos = 0.65*( np.array( [ points_x[1], points_y[1] ] ) - np.array( [ points_x[0], points_y[0] ] ) ) + np.array([0,0.17])
         plt.text(text_pos[0], text_pos[1], r'$b_1$', fontsize=12, color='blue')
      elif counter == 4 :
         plt.arrow(points_x[0], points_y[0], points_x[1], points_y[1], color='blue', linestyle='--', head_width=0.08, length_includes_head=True)
         
         text_pos = 0.75*( np.array( [ points_x[1], points_y[1] ] ) - np.array( [ points_x[0], points_y[0] ] ) ) + np.array([-0.22,0.0])
         plt.text(text_pos[0], text_pos[1], r'$b_2$', fontsize=12, color='blue')
      else :
         plt.plot(points_x, points_y, linewidth=1.5, color='black', linestyle='-', marker='', markersize=5 )

      points_x = []
      points_y = []
   else :
      points_x.append( float(data_str[0]) ) 
      points_y.append( float(data_str[1]) ) 
      # if counter == 1: print( points[-1] )

axes.set_xlim([ -box_lim, box_lim ])
axes.set_ylim([ -box_lim, box_lim ])

plt.legend()

outfile = "kpath_"+sys.argv[1]+".png"
print( 'kpath output file:', outfile )
plt.savefig(outfile)