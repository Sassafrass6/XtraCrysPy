from XtraCrysPy import XCP_BZ as XCP
import numpy as np

data = np.load('data_files/Fermi_surf_band_46_0.80x80x60.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))

# Two types of isosurfaces, one clipped and one unclipped

# Clipped
iso_vals_c = [-1.8, -1.2, -1, -.8]
colors_c = (250*np.array([[.6,.8,0], [0,.6,.3], [0,.5,1], [.5,0,.8]])).astype(int)

# Clipping planes
plane_centers = [[0,0,.05], [0,0,-.05]]
plane_normals = [[0,1,1], [0,-1,-1]]
clips = [plane_centers, plane_normals]

# Unclipped (with some transparency)
iso_vals_uc = [-2.4, -2.2]
colors_uc = (250*np.array([[1,0,0,.8], [1,.2,0,.7]])).astype(int)

xcp = XCP.XCP_BZ(model='data_files/Te_L.scf.in')
xcp.render_iso_surface(data, iso_vals=iso_vals_c, colors=colors_c, clip_planes=clips)
xcp.render_iso_surface(data, iso_vals=iso_vals_uc, colors=colors_uc, disp_all=True)
xcp.start_crystal_view()
