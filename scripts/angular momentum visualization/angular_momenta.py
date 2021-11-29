# Python 3 required due to Mayavi
# source python-select local3
# Prepare two inputcards for the spin-wave dispersion program
# One with a kpoint PATH and other for kpoint MESH 

from scipy.interpolate import griddata
import numpy as np
from mayavi import mlab

system = 'AF20x20c2x2spiral'; shift=0.0
system = 'AF20x20p2x1spiral'; shift=1.0
system = 'AF2x2c2x2'; shift=0.0

# num_modes = 40
filenameSpectrMESH =       "outputfiles/disp_unfol_"+system+"MESH.dat"
filenameAngMomMESH = "outputfiles/angular_momentum_"+system+"MESH.dat"
filenameAngMomPATH = "outputfiles/angular_momentum_"+system+"PATH.dat"

############################################
#     Reading data angular momenta PATH
############################################


file = open( filenameAngMomPATH, 'r' )
lines=file.readlines()

kxs = []
kys = []

mode_energies_allkpoints = []
mode_Svectors_allkpoints = []

for iline, line in enumerate(lines):
   data_str=line.split()

   #Skip the first line if empty
   if iline == 0 and len(data_str) == 0 : continue

   #Make sure there is something in the line
   if not len(data_str) == 0 : 

      #Read information of the kpoint
      if data_str[0] == 'kpointindex:' :

         kpointindex =  int(data_str[1])

         kxs.append( float(data_str[3]))
         kys.append( float(data_str[4]))

         #Prepare the array to store the energy and the Svectors for this kpoint
         mode_energies = []
         mode_Svectors = []

         # if kpointindex == 4 : break

         continue

      #For each kpoint, read information for every mode
      if data_str[0] == 'mode=' :
         mode_index = int(data_str[1])-1
         mode_energies.append(float(data_str[17]))
         mode_Svectors.append( [float(data_str[4]), float(data_str[5]), float(data_str[6])] )

   else : # if len(data_str) == 0 

      #After reading the info for every mode, store them
      mode_energies_allkpoints.append(mode_energies)
      mode_Svectors_allkpoints.append(mode_Svectors)

num_modes = mode_index + 1
num_kpoints = len(kxs)
print("number of K points in PATH:", num_kpoints)
print("number of modes:", num_modes)

# pp = mlab.points3d(kxs, kys, len(kxs)*[0], scale_factor=0.5, color=(0.0,0.5,1))
# pp = mlab.points3d([np.pi], [np.pi*0], [0], scale_factor=0.5, color=(0.0,0.5,1))

############################################
#     Reading data spectrum MESH
############################################


file = open( filenameSpectrMESH, 'r' )
lines=file.readlines()

omegas = []
spectrums = []
      
omegas_allkpoints = []
spectrums_allkpoints = []


for line in lines:
   data_str=line.split()

   # If the line is empty we store the omegas for a given kpoint
   if len(data_str) == 0 and len(omegas) > 0 : 
      omegas_allkpoints.append(omegas)
      spectrums_allkpoints.append(spectrums)

      omegas = []
      spectrums = []

      continue

   omegas.append( float(data_str[1]) )
   spectrum = np.array(data_str[2:])
   spectrum = spectrum.astype(np.float)

   spectrums.append(spectrum)


############################################
#    Plotting angular momenta PATH
############################################

energy_of_mode = np.zeros((num_kpoints,num_modes), dtype=np.float64)
Svectors_of_mode = np.zeros((num_kpoints,num_modes,3), dtype=np.float64)

for kpoint, (mode_energies,mode_Svectors) in enumerate(zip(mode_energies_allkpoints,mode_Svectors_allkpoints)) :
   for mode, (energy,Svector) in enumerate(zip(mode_energies,mode_Svectors)) :
      energy_of_mode[kpoint,mode] = energy
      Svectors_of_mode[kpoint,mode,:] = Svector


kxs = np.array(kxs)
kys = np.array(kys)
# print("position of the kpoint path")
# print(kxs)
# print(kys)

#The energy of all kpoints for mode
for mode in range(num_modes) :
   z = energy_of_mode[:,mode]
   Sx = Svectors_of_mode[:,mode,0]
   Sy = Svectors_of_mode[:,mode,1]
   Sz = Svectors_of_mode[:,mode,2]

   normS = []
   Sx_renormalized = []
   Sy_renormalized = []
   Sz_renormalized = []
   for xx,yy,zz in zip(Sx, Sy, Sz):
      normS.append(np.linalg.norm([xx,yy,zz]) )

      if normS[-1] >= 1.e-10 :
         Sx_renormalized.append(xx/normS[-1])
         Sy_renormalized.append(yy/normS[-1])
         Sz_renormalized.append(zz/normS[-1])
      else :
         Sx_renormalized.append(xx)
         Sy_renormalized.append(yy)
         Sz_renormalized.append(zz)


   Sx_renormalized = np.array(Sx_renormalized)
   Sy_renormalized = np.array(Sy_renormalized)
   Sz_renormalized = np.array(Sz_renormalized)

   theta = 0#np.pi/4.0
   #          -pi < arctan2 < pi
   # colors = (np.arctan2(Sx_renormalized*np.cos(theta)-Sy_renormalized*np.sin(theta),Sx_renormalized*np.sin(theta)+Sy_renormalized*np.cos(theta))/np.pi + 1.0)/2.0
   colors = ((np.arctan2(Sx_renormalized*np.cos(theta)-Sz_renormalized*np.sin(theta),Sx_renormalized*np.sin(theta)+Sz_renormalized*np.cos(theta))/np.pi + 1.0)/2.0)
   # print(colors)
   # exit()

   # p = mlab.points3d(kxs, kys, z, scale_factor=.2)
   s = mlab.plot3d(kxs, kys, z,normS, tube_radius=0.02, color=(0,0,0), transparent=True, opacity=0.1)# ,normS, colormap='YlOrBr')
   obj = mlab.quiver3d(kxs -0* Sx/4, kys -0* Sy/4, z -0* Sz/4, Sx, Sy, Sz, line_width=10, scale_factor=0.5, mode='arrow',resolution=25) #colormap='gist_rainbow'
   obj.glyph.color_mode = "color_by_scalar"
   obj.mlab_source.dataset.point_data.scalars = colors

############################################
#     Reading data kpoint MESH
############################################

file = open( filenameAngMomMESH, 'r' )
lines=file.readlines()

kxs = []
kys = []

for line in lines:
   data_str=line.split()

   #Skip the loop if the line is empty
   if len(data_str) == 0 : continue

   #Read information of the kpoint
   if data_str[0] == 'kpointindex:' :

      kxs.append( float(data_str[3]))
      kys.append( float(data_str[4]))

      continue

num_kpoints = len(kxs)
print("number of K points in MESH:", num_kpoints)

############################################
#    Plotting spectrum 3D
############################################

num_omega = len(omegas_allkpoints[0])
print('num_omega:', num_omega)

omega_of_mode    = np.zeros((num_kpoints,num_omega)   , dtype=np.float64)
spectrum_of_mode = np.zeros((num_kpoints,num_omega,12), dtype=np.float64)

for kpoint, (omegas,spectrums) in enumerate(zip(omegas_allkpoints,spectrums_allkpoints)) :
   for omega_index, (omegas,spectrum) in enumerate(zip(omegas,spectrums)) :
      omega_of_mode[kpoint,omega_index] = omegas
      spectrum_of_mode[kpoint,omega_index,:] = spectrum


dimension = int(np.sqrt(len(kxs)))
print('kx mesh:', dimension, 'x', dimension )
kxs2ds = np.reshape(kxs, (-1, dimension))
kys2ds = np.reshape(kys, (-1, dimension))
kxsgrid = []
kysgrid = kys2ds[0]+np.pi*shift
for kxs2d in kxs2ds :
   kxsgrid.append(kxs2d[0]-np.pi*shift)

kxsgrid = np.array(kxsgrid)
kysgrid = np.array(kysgrid)

x = kxsgrid
y = kysgrid
z = omega_of_mode[0,:]

# print("position of the kpoint mesh")
# print(x)
# print(y)
# print(z)
print('len x, y, z:', len(x), len(y), len(z))
# exit()
y, x, z = np.meshgrid(x, y, z)

omegas_ky = []
kpoint = -1
for j, kx in enumerate(kxsgrid) :
   omegas_kx = []
   for ky in kysgrid :
      omegas = []
      kpoint += 1
      #For every omega
      for oidx, (omega, spectrum_array) in enumerate(zip(omega_of_mode[kpoint,:], spectrum_of_mode[kpoint,:])) :  
         # print(j, kpoint, oidx,'kx',kx,'ky',ky, 'omega', omega, 'spec', spectrum_array[0])

         omegas.append(sum(spectrum_array[0:4]))
         # omegas.append(spectrum_array[0])

      omegas_kx.append(omegas)

   omegas_ky.append(omegas_kx)

maxvalue = np.amax(omegas_ky)
print("max. spectrum value:", maxvalue)
scalars = omegas_ky


# obj = mlab.contour3d(x,y,z,scalars,contours=20, transparent=True, vmin=maxvalue*0., vmax=maxvalue*1.)

# vol=mlab.pipeline.volume(mlab.pipeline.scalar_field(x,y,z,scalars), vmin=maxvalue*0.0, vmax=maxvalue*0.9)

# mlab.volume_slice(x,y,z,scalars, vmin=maxvalue*0.2, vmax=maxvalue*0.9)

# slide1=mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(x,y,z,scalars),
#                             plane_orientation='x_axes',
#                             slice_index=dimension/2-1,
#                             vmin=0.0, vmax=maxvalue*0.9,
#                         )

# slide2=mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(x,y,z,scalars),
#                             plane_orientation='y_axes',
#                             slice_index=dimension/2-1,
#                             vmin=0.0, vmax=maxvalue*0.9,
#                         )

# slide1.edit_traits()
# slide2.edit_traits()

# mlab.show_pipeline()
mlab.show()













