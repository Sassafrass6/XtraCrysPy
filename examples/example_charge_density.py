from XtraCrysPy.Atomic import Atomic
from XtraCrysPy.file_io import read_XSF

if '__main__' == __name__:

  a_info = {'colors':{'Si':(0,.2,1)},
            'bonds':{'distance':3.3, 'thickness':0.5}}

  fname = 'data_files/Si.scf.in'
  xcp = Atomic(model=fname, params=a_info, nsc=(3,3,2))

  density_file = 'data_files/Si.density_0.xsf'
  data = read_XSF(density_file)

  xcp.render_iso_surface(data[0], iso_vals=.000075, colors=[.9,.2,0,.8])

  # The camera's position and orientation is set manually.
  #  These coordinates can be obtained any time by pressing 'o'+SHIFT
  cam_pos = [56.0194, -56.0907, 18.4797]
  cam_foc = [-12.9433, 12.9455, 15.5912]
  cam_up = [0.5538, 0.5783, 0.5990]

  xcp.start_crystal_view(camera_pos=cam_pos,
                         camera_focal=cam_foc,
                         camera_up=cam_up)

