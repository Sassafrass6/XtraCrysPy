from XtraCrysPy import XCP_BZ as XCP
import numpy as np

data = np.load('data_files/Al_Fermi_surf_band_2_0.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))

iso_vals = 4
colors = [1,0,0]
xcp = XCP.XCP_BZ(model='data_files/Al.scf.in')
xcp.render_iso_surface(data, iso_vals=iso_vals, colors=colors)
xcp.start_crystal_view()
