from XtraCrysPy import XCP_BZ as XCP
import numpy as np

# Te-L scf file
xcp = XCP.XCP_BZ(model='./data_files/Te_L.scf.in', background=(0,0,0))

# Load Fermi surface and center grid with fftshift
data = np.load('./data_files/Te_L.Fermi_surf_band_47_0.npz')['nameband']
data = np.fft.fftshift(data, axes=(0,1,2))

# Load spin texture and center grid with fftshift
spins = np.load('./data_files/Te_L.spin_text_band_47_0.npz')['spinband']
spins = np.fft.fftshift(np.real(spins), axes=(0,1,2))

# Had to roll the spin data by -1 on each axis to get the correct centering.
#   Unsure if this effect is specific to the data or plotting routines...
spins = np.roll(np.roll(np.roll(spins,-1,axis=0),-1,axis=1),-1,axis=2)

# Select iso-value and colors
#   Surface can have opacity moduleated (RGB or RGBA)
#   Arrows only support 3-colors currently (RGB)
iso = -.035
arrow_color = np.array([255,125,0])
surface_color = np.array([0,125,255,175])
trans_color = np.array([0,0,0,0])

# "Inflate" the BZ clipping planes to allow plotting of the pockets slightly outside of the BZ
planes = 1.2*xcp.bound_planes.copy()

# Make cut planes close to the BZ boundary
z_bound_top = xcp.bound_planes.copy()[0,0]
dx = 1e-5*np.array([0,0,1])
# Cut centers and normals
cuts = np.array([[z_bound_top+dx, z_bound_top-dx], [[0,0,1],[0,0,-1]]])

# Disable the default boundary clipping and provide the updated planes
xcp.render_iso_surface(data,
                       iso_vals=iso,
                       colors=surface_color,
                       clip_boundary=False,
                       clip_planes=planes)

xcp.render_iso_surface(data,
                       iso_vals=iso,
                       colors=surface_color,
                       arrows=spins,
                       arrow_scale=0.01,
                       arrow_colors=arrow_color,
                       clip_boundary=False,
                       clip_planes=cuts,
                       disp_all=True)

xcp.start_crystal_view()
