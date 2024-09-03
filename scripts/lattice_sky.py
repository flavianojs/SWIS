#!/usr/bin/env python

import numpy as np
import gr, gr3
import colorsys
import math
import sys

########################################################################
#                 file with vector field data
########################################################################
spirit = False

lattice_shift_x = 0 #10.5
lattice_shift_y = 0 #10.5

# Get filename from command line
filename = sys.argv[1]
file = open(filename)
lines_file_spinconfig=file.readlines()

if spirit :
    # file2   = open("inputfiles/lattice8x8.dat")
    filename = sys.argv[2]
    file2   = open(filename)
    lines_file_lattice = file2.readlines()

filename = sys.argv[3]



# Set some arrays
# Coordinates of the vectors
x=[]
y=[]
z=[]
# Vectors
phi=[]
theta=[]
nr=[]

# line = lines_file_spinconfig[0]
# if line.strip().startswith('#') :
#     spirit = True

# Read data columns into the respective arrays
if not spirit :
    for line in lines_file_spinconfig[0:]:
        if not line.strip().startswith('#') :
            data_str=line.split()
            x.append(float(data_str[0]))
            y.append(float(data_str[1]))
            z.append(float(data_str[2]))
            phi.append(float(data_str[3]))
            theta.append(float(data_str[4]))
            nr.append(float(data_str[5]))
else :
    for line in lines_file_lattice[0:]:
        if not line.strip().startswith('#') :
            data_str=line.split()
            x.append(float(data_str[0])-lattice_shift_x)
            y.append(float(data_str[1])-lattice_shift_y)
            z.append(float(data_str[2]))
    
    for line in lines_file_spinconfig[0:]:
        if not line.strip().startswith('#') :
            data_str=line.split()
            phi.append(float(data_str[0]))
            theta.append(float(data_str[1]))
            nr.append(float(data_str[2]))



positions = np.zeros((len(x), 3), dtype=np.float32)
positions[:, 0] = x
positions[:, 1] = y
positions[:, 2] = z

directions = np.zeros(positions.shape, dtype=np.float32)
if spirit :
    directions[:, 0] = phi
    directions[:, 1] = theta
    directions[:, 2] = nr
else :
    directions[:, 0] = np.sin(theta)*np.cos(phi)*nr
    directions[:, 1] = np.sin(theta)*np.sin(phi)*nr
    directions[:, 2] = np.cos(theta)*nr
    # directions[:, 2] = np.sqrt(theta)

colors = np.zeros(positions.shape, dtype=np.float32)

for i, direction in enumerate(directions):
    x = np.sqrt(np.square(direction[0])+np.square(direction[1]))
    hue =   np.arctan2(x,direction[2]) / np.pi / 1.5
    saturation = 1
    value = 1

    colors[i] = colorsys.hsv_to_rgb(hue, saturation, value)

rescale = 0.5

aa = 0.6
bb = 0.7
for i, position in enumerate(positions):
    # if i < 64 : 
    if -aa < position[0] < aa and -bb < position[1] < bb: # for the structure vs B field
        # print (position[0])
        size_arrow = np.linalg.norm( directions[i])
        gr3.drawspins(4*positions[i], directions[i], colors[i],
                      cone_radius=0.28*rescale*size_arrow, cylinder_radius=0.15*rescale*size_arrow,
                      cone_height=0.6*rescale*size_arrow, cylinder_height=0.7*rescale*size_arrow)


# #To plot a axis system
# positions = np.zeros((3, 3), dtype=np.float32)
# positions[0,:] = [15, 0, 0]
# positions[1,:] = [15, 0, 0]
# positions[2,:] = [15, 0, 0]
# axis = np.zeros(positions.shape, dtype=np.float32)
# axis[0,:] = [ 5, 0, 0]
# axis[1,:] = [ 0,15, 0]
# axis[2,:] = [ 0, 0,45]

# colors = np.zeros(positions.shape, dtype=np.float32)
# for i, axi in enumerate(axis):
#     x = np.sqrt(np.square(axi[0])+np.square(axi[1]))
#     hue =   np.arctan2(x,axi[2]) / np.pi / 1.5
#     saturation = 1
#     value = 1
#     colors[i] = colorsys.hsv_to_rgb(hue, saturation, value)

# gr3.drawspins(positions, axis, colors,
#               cone_radius=0.25, cylinder_radius=0.15,
#               cone_height=0.6, cylinder_height=0.4)


# gr3.setprojectiontype(1)
# gr.setperspectiveprojection(30, -30, 0)
# gr3.setorthographicprojection(float left, float right, float bottom, float top, float znear, float zfar))
# gr3.setorthographicprojection(-30, 30, -30, 30, 30, -30)

# gr3.setbackgroundcolor(0.2, 0.2, 0.2, 1.0)
# gr3.setbackgroundcolor(0., 0., 0., 1.)
gr3.setbackgroundcolor(1. , 1. , 1. , 0.5)

# gr3.cameralookat(0,-26, 26, 0, 0, 0, 0, 1, 0) #Angle 2: 
gr3.cameralookat(0, 0, 30, 0, 0, 0, 0, 1, 0) #Angle 1: Top view
gr3.cameralookat(0, 0, 30, 0, 0, 0, 0, 1, 0) #Angle 1: Top view

gr3.setlightdirection(1,0,2)
gr3.setlightdirection(0,0,1)
gr3.export(filename, 2000, 2000)
# gr3.export(filename, 1500, 2000) # for the structure vs B field

