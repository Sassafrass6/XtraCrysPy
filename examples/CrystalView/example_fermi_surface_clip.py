from XtraCrysPy import XCP_BZ as XCP
import numpy as np

data = np.load('data_files/Fermi_surf_band_46_0.80x80x60.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))

iso_vals = np.linspace(-2.4, -.5, 5)
colors = (250*np.array([[1,0,0,1], [1,.4,0,.8], [.6,.8,0,.7], [0,1,.2,.6], [0,.5,1,.4]])).astype(int)

plane_centers = [[0,0,.02], [0,0,-.02]]
plane_normals = [[0,1,1], [0,-1,-1]]
clips = [plane_centers, plane_normals]

xcp = XCP.XCP_BZ(model='data_files/Te_L.scf.in')
xcp.render_iso_surface(data, iso_vals=iso_vals, colors=colors, disp_all=True, clip_planes=clips)
xcp.start_crystal_view()
