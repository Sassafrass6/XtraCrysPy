from XtraCrysPy.Atomic import Atomic
import numpy as np

N = 36
# Generate some positions for each of two particles
positions = np.empty((N,2,3), dtype=float)

domain = np.linspace(0, 2*np.pi, N)
for i,theta in enumerate(domain):
  x,y = np.cos(theta), np.sin(theta)
  positions[i,0,:] = [x, y, 0]
  positions[i,1,:] = [-x, -y, 0]

# Color one red and the other blue
params = {'abc':positions,
          'species':['s1', 's2'],
          'colors':{'s1':(1,0,0),
                    's2':(0,0,1)}}

# Render the multi-frame visualization
xcp = Atomic(params=params, multi_frame=True)

# Animate with frame duration 30ms, 1 step per frame, and a restart delay of 20 frames.
xcp.animate(fdt=30, spf=1, restart_delay=20)
xcp.start_crystal_view()
