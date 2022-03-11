from XtraCrysPy import XCP_BZ as XCP
import numpy as np


def read_bxsf ( fname:str ):

  lines = None
  with open(fname, 'r') as f:
    lines = f.readlines()

  ng = np.array([int(v) for v in lines[10].replace('\n','').split()])
  rlat = np.empty((3,3), dtype=float)
  for i in range(3):
    rlat[i,:] = np.array([float(v) for v in lines[12+i].replace('\n','').split()])

  ind = 15
  ng[0] = 2
  data = np.empty(ng, dtype=float)
  for i in range(ng[0]):
    ind += 1
    for j in range(ng[1]):
      for k,v in enumerate(lines[ind].split()):
        data[i,j,k] = float(v)
        ind += 1
  return rlat,data
  
#eigs = np.load('data_files/SnTe.fermi_surf_band_4.npz')['nameband']
#rlat,data = read_bxsf('data_files/SnTe.spin_berry_z_xy.bxsf')
data = np.load('data_files/Fermi_surf_band_46_0.80x80x60.npz')['nameband']

colors = np.zeros((*data.shape, 3))
for i in range(data.shape[2]):
  colors[:,:,i,:] = np.array([150,0,255*i/data.shape[2]])

xcp = XCP.XCP_BZ(model='data_files/scf-last.in')
xcp.render_iso_surface(data, iso_val=-1, colors=colors)
xcp.start_crystal_view()
