from XtraCrysPy.Atomic import Atomic
import numpy as np

# Generate a radially increasing data grid.
# Color Green-Blue from (-1,0) on Z, and Red-Green on (0,1).
# Determine arrow directions with circulation about Z.
N = 64
grid = np.empty((N,N,N), dtype=float)
col = np.empty((N,N,N,3), dtype=float)
arrows = np.empty((N,N,N,3), dtype=float)
dg = np.linspace(-1,1,N)
for i,x in enumerate(dg):
  for j,y in enumerate(dg):
    for k,z in enumerate(dg):
      col[i,k,j,:] = [0,1,0]
      col[i,k,j,(2 if z<0 else 0)] = np.abs(z)
      grid[i,j,k] = np.linalg.norm([x,y,2*z])
      if x == 0 and y == 0:
        arrows[i,j,k,:] = [0,0,1]
      elif x == 0:
        arrows[i,j,k,:] = [(1 if y>0 else -1), 0, 0]
      else:
        th = np.arctan(y/x)
        if x < 0:
          th += np.pi
        arrows[i,j,k,:] = [np.sin(th), -np.cos(th), 0]

# Create XCP object and render iso-surface with texture.
xcp = Atomic(boundary=False, background=(1,1,1))
xcp.render_iso_surface(grid, iso_vals=0.7, colors=col, arrows=arrows, arrow_scale=2, arrow_colors=[0,0,1], arrow_anchor='mid', arrow_spacing=0.033, clip_boundary=True)
xcp.start_crystal_view()

