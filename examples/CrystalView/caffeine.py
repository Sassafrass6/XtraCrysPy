from XtraCrysPy import XtraCrysPy as XCP
import numpy as np

if '__main__' == __name__:

  # Define params
  a = 20.

  # Unit origin
  origin = [0,0,0]

  # Cell parameters (Normalized)
  fcc_vecs = [[a,0,0], [0,a,0], [0,0,a]]

  # Atomic basis
  atoms = 2*(np.array([[0.47,-3.1271,-0.9686,2.2182,-1.3477,1.4119,0.8579,0.3897,0.0307,-1.9061,2.5032,-1.4276,3.1926,-2.2969,3.5163,-1.0451,-2.5186,-1.0447,4.1992,3.0468,3.0466,-1.8087,-2.9322,-2.9346],[2.5688,-0.4436,-1.3125,0.1412,1.0797,-1.9372,0.2592,-1.0264,1.422,-0.2495,-1.1998,-2.696,1.2061,2.1881,-1.5787,-3.1973,-2.7596,-3.1963,0.7801,1.8092,1.8083,3.1651,2.1027,2.1021],[0.0006,-0.0003,0,-0.0003,-0.0001,0.0002,-0.0008,-0.0004,-0.0006,-0.0004,0.0003,0.0008,0.0003,0.0007,0.0008,-0.8937,0.0011,0.8957,0.0002,-0.8992,0.9004,-0.0003,0.8881,-0.8849]]).T+5)

  # Label for each atom and colors
  labels = ['O']*2 + ['N']*4 + ['C']*8 + ['H']*10

  # Bond distance between Ga and As, in ????
  bonds = {'C_C':1.46*2, 'C_N':1.46*2, 'C_O':1.47*2, 'C_H':1.1*2}

  # Species information, including radius and color

  species = {'O':{'color':(1,0,0,1),'radius':.7}, 'C':{'color':(0.1,0.1,0.1,1),'radius':.8}, 'N':{'color':(0,0,1,1),'radius':.8},'H':{'color':(1,1,1,1),'radius':.5}}

  
  # Construct model with XtraCrysPy
  xcpy = XCP.XtraCrysPy(lattice=fcc_vecs, basis=atoms, basis_labels=labels, species=species, bonds=bonds, origin=origin)

  # Write an XML file for blender
  #xcpy.write_blender_xml(fname='Caffeine_Ball.xml')

  xcpy.start_cryspy_view()
