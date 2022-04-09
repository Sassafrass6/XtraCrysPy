from XtraCrysPy.Atomic import Atomic
from XtraCrysPy.file_io import read_XSF

if '__main__' == __name__:

  a_info = {'colors':{'Si':(0,.2,1)},
            'bonds':{'distance':6.5, 'thickness':0.5}}

  fname = 'data_files/Si.scf.in'
  xcp = Atomic(model=fname, params=a_info, nsc=(3,3,2))

  density_file = 'data_files/Si.density_0.xsf'
  data = read_XSF(density_file)

  xcp.render_iso_surface(data[0], iso_vals=.000075, colors=[.9,.2,0,.8])
  xcp.start_crystal_view()

