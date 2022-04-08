from XtraCrysPy import Atomic,Model
import numpy as np

# Define atomic parameters explicitly
struct = {'lattice':5.4*np.array([[-1,0,1],[0,1,1],[-1,1,0]]),
          'species':['Si','Ge'],
          'abc':np.array([[0.,0,0],[.25,.25,.25]]),
          'colors':{'Si':[1,0,0], 'Ge':[0,0,1]},
          'radii':{'Si':1.2, 'Ge':1.2},
          'bonds':{'Si_Ge':6}}

xcp = Atomic.Atomic(params=struct, nsc=(2,2,2), bond_type='Primary', background=(.5,.5,.5), sel_type='Distance')
xcp.start_crystal_view()

