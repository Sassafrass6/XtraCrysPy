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
#data = np.load('data_files/Fermi_surf_band_46_0.80x80x60.npz')['nameband']
data = np.load('data_files/Fermi_surf_band_103_0.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))
spins = np.load('data_files/spin_text_band_103_0.npz')['spinband']
spins = np.fft.fftshift(spins, axes=(0,1,2))
spins = np.real(spins)

data = data[::2,::2,::2]
spins = spins[::2,::2,::2]
ds = data.shape
colors = np.zeros((*ds, 4))
for i in range(ds[0]):
  for j in range(ds[2]):
    a = i if i<ds[0]//2 else ds[0]-i
    b = j if j<ds[2]//2 else ds[2]-j
    col = [a/ds[0], b/ds[2], .5, .9]
    colors[i,:,j,:] = 255 * np.array(col)

plane_centers = [[0,0,.01], [0,0,-.01]]
plane_normals = [[0,.4,1], [0,-.4,-1]]
clips = [plane_centers, plane_normals]
iso_vals = np.linspace(-1.4, -.2, 6)
xcp = XCP.XCP_BZ(model='data_files/Te_L.scf.in')
xcp.render_iso_surface(data, iso_vals=iso_vals[::2], colors=colors, arrows=spins, arrow_colors=(1,0,0), clip_planes=clips)
xcp.render_iso_surface(data, iso_vals=iso_vals[1::2], colors=colors, arrows=None, clip_planes=clips, disp_all=True)
xcp.start_crystal_view()
