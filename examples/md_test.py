from XtraCrysPy.MolDyn import MolDyn
import numpy as np

positions = np.empty((3,3), dtype=float)
positions[0] = np.array([-1,0,0])
positions[1] = np.array([0,0,0])
positions[2] = np.array([.8,.8,.8])

# Color one red and the other blue
params = {'abc':positions,
          'species':['s1', 's2', 's3'],
          'radii':{'s1':0.5,'s2':0.5,'s3':0.5},
          'colors':{'s1':(1,0,0),
                    's2':(0,1,0),
                    's3':(0,0,1)}}

# Render the multi-frame visualization
xcp = MolDyn(params=params)
xcp.start_crystal_view()
