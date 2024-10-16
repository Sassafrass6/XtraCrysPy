from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  # 2D SnTe with colors red & blue for Sn & Te respectively.
  # Bonds are drawn between Sn & Te atoms if they are closer together than 6 a.u.
  a_info = {'colors':{'Sn':(1,0,0),'Te':(0,0,1)}, 'radii':{'Sn':1,'Te':1.1}, 'bonds':{'Sn_Te':3.3}}
  fname = 'data_files/SnTe.scf.in'
  xcp = Atomic(model=fname, params=a_info, boundary=False, nsc=(3,3,1))
  xcp.start_crystal_view()

