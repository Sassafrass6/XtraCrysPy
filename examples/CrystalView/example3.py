from XtraCrysPy import XCP_Atoms,Model

if '__main__' == __name__:

  # 2D SnTe with colors red & blue for Sn & Te respectively.
  # Bonds are drawn between Sn & Te atoms if they are closer together than 6 angstrom
  a_info = {'colors':{'Sn':(1,0,0),'Te':(0,0,1)}, 'radii':{'Sn':1,'Te':1.1}}
  fname = 'SnTe_data/SnTe.scf.in'
  xcp = XCP_Atoms.XCP_Atoms(model=fname, params=a_info, axes=False, nsc=(3,3,1))
  xcp.start_crystal_view()

