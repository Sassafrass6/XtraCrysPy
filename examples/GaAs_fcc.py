from XtraCrysPy import Atomic,Model
import numpy as np

# Define atomic parameters explicitly
struct = {'species':['Ga','As'],
          'lunit':'bohr',
          'lattice':5.4*np.array([[-1,0,1],[0,1,1],[-1,1,0]]),
          'aunit':'crystal',
          'abc':np.array([[0.,0,0],[.25,.25,.25]]),
          'colors':{'Ga':[1,0,0], 'As':[0,0,1]},
          'radii':{'Ga':1.2, 'As':1.2},
          'bonds':{'Ga_As':6}}

xcp = Atomic.Atomic(params=struct, nsc=(2,2,2), bond_type='Primary', background=(.5,.5,.5), sel_type='Distance')
xcp.start_crystal_view()

