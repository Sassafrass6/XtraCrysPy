from XtraCrysPy import XCP_BZ as XCP
import numpy as np

data = np.load('data_files/Fermi_surf_band_103_0.npz')['nameband']
data = data[::2,::2,::2]
data = np.fft.fftshift(data, axes=(0,1,2))

spins = np.load('data_files/spin_text_band_103_0.npz')['spinband']
spins = spins[::2,::2,::2]
spins = np.fft.fftshift(spins, axes=(0,1,2))
spins = np.real(spins)

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

xcp = XCP.XCP_BZ(model='data_files/Te_L.scf.in', background=(1,1,1))
xcp.render_iso_surface(data, iso_vals=iso_vals[::2], colors=colors, arrows=spins, arrow_colors=(1,0,0), clip_planes=clips)
xcp.render_iso_surface(data, iso_vals=iso_vals[1::2], colors=colors, arrows=None, clip_planes=clips, disp_all=True)
xcp.start_crystal_view()
