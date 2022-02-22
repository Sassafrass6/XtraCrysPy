from XtraCrysPy import XtraCrysPy as XCP
import numpy as np

alat = 4
a_vec = alat * np.array([[1,0,0],[np.sqrt(3)/2,.5,0],[0,0,1.3]])

omega = a_vec[0] @ np.cross(a_vec[1], a_vec[2])
print(omega)
b_vec = np.zeros_like(a_vec)
b_vec[0,:] = 2*np.pi/omega * np.cross(a_vec[1], a_vec[2])
b_vec[1,:] = 2*np.pi/omega * np.cross(a_vec[2], a_vec[0])
b_vec[2,:] = 2*np.pi/omega * np.cross(a_vec[0], a_vec[1])

print(b_vec)

xcp = XCP.XtraCrysPy()
xcp.start_cryspy_view()
xcp.view.draw_BZ_boundary(b_vec=b_vec)

