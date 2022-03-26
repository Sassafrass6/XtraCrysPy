from XtraCrysPy import XCP_BZ as XCP
import numpy as np

data = np.load('data_files/Fermi_surf_band_46_0.80x80x60.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))

ds = data.shape
colors = np.zeros((*ds, 3))
for i in range(ds[0]):
  for j in range(ds[2]):
    a = i if i<ds[0]//2 else ds[0]-i
    b = j if j<ds[2]//2 else ds[2]-j
    col = [a/ds[0], b/ds[2], .5]
    colors[i,:,j,:] = 255 * np.array(col)

iso_vals = np.linspace(-1.3, 0, 5)
xcp = XCP.XCP_BZ(model='data_files/Te_L.scf.in')
xcp.render_iso_surface(data, iso_vals=iso_vals, colors=colors)
xcp.start_crystal_view()
