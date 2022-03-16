#!/usr/bin/env python

import numpy as np
import math

#
# Input section
#

toprint = False

#Interaction file to be converted
filename   = "inputfiles/pair_MnSi_rad4.5.txt"
outputfile = "inputfiles/pair_MnSi_rad4.5rect.txt"

#Initial basis
num_atom_basis_ini = 8

a1_ini       = np.array([ 1.0             ,  0.0             ,  0.0 ])
a2_ini       = np.array([-1.0             ,  1.73205080756888,  0.0 ])
a3_ini       = np.array([ 0.0             ,  0.0             ,  0.69679557694152 ])

basis_ini    = np.zeros((num_atom_basis_ini,3),      dtype=np.float64)
basis_ini[0] = np.array([ 0.000000,   0.382000,   0.750000 ])
basis_ini[1] = np.array([ 0.764000,   0.382000,   0.250000 ])
basis_ini[2] = np.array([ 0.000000,   0.118000,   0.250000 ])
basis_ini[3] = np.array([ 0.236000,   0.118000,   0.750000 ])
basis_ini[4] = np.array([ 0.764000,   0.882000,   0.250000 ])
basis_ini[5] = np.array([ 0.000000,   0.882000,   0.750000 ])
basis_ini[6] = np.array([ 0.236000,   0.618000,   0.750000 ])
basis_ini[7] = np.array([ 0.000000,   0.618000,   0.250000 ])

#Final basis
num_atom_basis_fin = 8

a1_fin       = np.array([ 1.0             ,  0.0             ,  0.0 ])
a2_fin       = np.array([ 0.0             ,  1.73205080756888,  0.0 ])
a3_fin       = np.array([ 0.0             ,  0.0             ,  0.69679557694152 ])

basis_fin    = np.zeros((num_atom_basis_fin,3),      dtype=np.float64)
basis_fin[0] = np.array([  -0.382000+0.5,   0.3820000+0.5,   0.750000 ])
basis_fin[1] = np.array([   0.382000+0.5,   0.3820000+0.5,   0.250000 ])
basis_fin[2] = np.array([  -0.118000+0.5,   0.1180000+0.5,   0.250000 ])
basis_fin[3] = np.array([   0.118000+0.5,   0.1180000+0.5,   0.750000 ])
basis_fin[4] = np.array([  -0.118000+0.5,  -0.1180000+0.5,   0.250000 ])
basis_fin[5] = np.array([   0.118000+0.5,  -0.1180000+0.5,   0.750000 ])
basis_fin[6] = np.array([  -0.382000+0.5,  -0.3820000+0.5,   0.750000 ])
basis_fin[7] = np.array([   0.382000+0.5,  -0.3820000+0.5,   0.250000 ])

# a1_fin       = np.array([ 1.0             ,  0.0             ,  0.0 ])
# a2_fin       = np.array([-1.0             ,  1.73205080756888,  0.0 ])
# a3_fin       = np.array([ 0.0             ,  0.0             ,  0.69679557694152 ])

# basis_fin    = np.zeros((num_atom_basis_fin,3),      dtype=np.float64)
# basis_fin[0] = np.array([ 0.000000,   0.382000,   0.750000 ])
# basis_fin[1] = np.array([ 0.764000,   0.382000,   0.250000 ])
# basis_fin[2] = np.array([ 0.000000,   0.118000,   0.250000 ])
# basis_fin[3] = np.array([ 0.236000,   0.118000,   0.750000 ])
# basis_fin[4] = np.array([ 0.764000,   0.882000,   0.250000 ])
# basis_fin[5] = np.array([ 0.000000,   0.882000,   0.750000 ])
# basis_fin[6] = np.array([ 0.236000,   0.618000,   0.750000 ])
# basis_fin[7] = np.array([ 0.000000,   0.618000,   0.250000 ])


print basis_fin
mapping = [[0], [1], [2], [3], [4], [5], [6], [7]]

#
# End Input section
#


file = open(filename)
lines=file.readlines()

target = open(outputfile, 'w')
target.write( lines[0] )  

# Set some arrays
# Coordinates of the vectors
# xs = []
# ys = []
# zs = []
ii = []
jj = []
das= []
dbs= []
dcs= []
# Vectors
Dxs = []
Dys = []
Dzs = []
DD = []
JJs = []

# RR = []
# positions = []

counts = np.zeros(num_atom_basis_ini, dtype=np.float64)

# Read data columns into the respective arrays
for line in lines[0:]:
   data_str=line.split()
   if data_str[0]!="i" :
      ii.append(  int(data_str[0 ]) )
      jj.append(  int(data_str[1 ]) )
      das.append( int(data_str[2 ]) )
      dbs.append( int(data_str[3 ]) )
      dcs.append( int(data_str[4 ]) )
      counts[ii[-1]] += 1

      JJs.append( float(data_str[5]) )
      DD.append(  float(data_str[6]) )
      Dxs.append( float(data_str[7]) )
      Dys.append( float(data_str[8]) )
      Dzs.append( float(data_str[9]) )

print 'Initial number of pair interaction in file:', len(DD)
print 'Number of pair for each atom in the initial unit cell:', ' %i '*num_atom_basis_ini%tuple(counts)

def cartesian_basis(basis,a1,a2,a3) :

   cartesian_basis = basis[0]*a1 + basis[1]*a2 + basis[2]*a3
   return cartesian_basis

num_int_pair_fin = 0
for atom in range(num_atom_basis_fin) :
   
   atom_basis = cartesian_basis(basis_fin[atom], a1_fin,a2_fin,a3_fin) 
   # print 'Calculating for atom ', atom, ' of the new basis.'

   num_interaction = 0
   # Loop over pair file to identify which pairs fits to the atoms in the new basis 
   # for i, j, da, db, dc, position, Dx, Dy, Dz, JJ, x, y, z in zip(ii, jj, das, dbs, dcs, positions, Dxs, Dys, Dzs, JJs, xs, ys, zs) :
   for i, j, da, db, dc, Dx, Dy, Dz, JJ, D in zip(ii, jj, das, dbs, dcs, Dxs, Dys, Dzs, JJs, DD) :
      if i == mapping[atom][0] :
         # if atom == 1 and i==2 and j==3 : toprint = True
      
         R_j  = da*a1_ini + db*a2_ini + dc*a3_ini + cartesian_basis(basis_ini[j], a1_ini,a2_ini,a3_ini)  
         R_ij = R_j - cartesian_basis(basis_ini[i], a1_ini,a2_ini,a3_ini) + atom_basis

         if toprint :
            print '%4i'*5%( i, j, da, db, dc), '%8.4f'*3%tuple(R_ij)

         da_aux = np.dot( a3_fin, np.cross(a2_fin, R_ij) ) / np.dot( a3_fin, np.cross(a2_fin, a1_fin) )
         db_aux = np.dot( a3_fin, np.cross(a1_fin, R_ij) ) / np.dot( a3_fin, np.cross(a1_fin, a2_fin) )
         dc_aux = np.dot( a1_fin, np.cross(a2_fin, R_ij) ) / np.dot( a1_fin, np.cross(a2_fin, a3_fin) )

         da_fin = int(np.floor( da_aux ))
         db_fin = int(np.floor( db_aux ))
         dc_fin = int(np.floor( dc_aux ))
         
         if( abs(da_aux-round(da_aux)) < 0.01 ): da_fin = int(round( da_aux ))
         if( abs(db_aux-round(db_aux)) < 0.01 ): db_fin = int(round( db_aux ))
         if( abs(dc_aux-round(dc_aux)) < 0.01 ): dc_fin = int(round( dc_aux ))

         found = False
         for atom2 in range(num_atom_basis_fin) :
            if j in mapping[atom2] :

               R_jaux = da_fin*a1_fin + db_fin*a2_fin + dc_fin*a3_fin + cartesian_basis(basis_fin[atom2], a1_fin,a2_fin,a3_fin)
               R_ijaux = R_jaux
               
               if toprint :
                  print '   - ', ' -'+'%4i'*3%( da_fin, db_fin, dc_fin), '%8.4f'*3%tuple(R_ijaux)

               diff = np.linalg.norm(R_ij-R_ijaux)
               if diff < 0.00001 :
                  i_fin = atom 
                  j_fin = atom2
                  found = True
                  R_ijaux2 = R_ijaux
                  
                  target.write( "%4i"*2%( i_fin, j_fin ) + "%4i"*3%( da_fin, db_fin, dc_fin ) + "%16.8f"*5%(JJ, D, Dx, Dy, Dz,  ) + "\n"  )
                  num_interaction += 1

         if found :
            if toprint :
               print '_______________________________________________'
               print '%4i'*5%( i, j, da, db, dc), '%8.4f'*3%tuple(R_ij)
               print '%4i'*5%( i_fin, j_fin, da_fin, db_fin, dc_fin), '%8.4f'*3%tuple(R_ijaux2)

         else :
            if toprint :
               print '_______________________________________________'
               print "didn't find"

         # if atom == 1 and i==2 and j==3 :
         if toprint : wait = input("type something: ")

   print 'Number of interaction for atom', atom, ':', num_interaction
   num_int_pair_fin += num_interaction

print 'Final number of pair interaction in file:', num_int_pair_fin
print 'Does it compare to the initial number???:', num_int_pair_fin * num_atom_basis_ini / num_atom_basis_fin












