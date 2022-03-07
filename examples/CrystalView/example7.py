from XtraCrysPy import XCP_Atoms,Model

if '__main__' == __name__:

  # 2D SnTe with colors red & blue for Sn & Te respectively.
  # Bonds are drawn between Sn & Te atoms if they are closer together than 6 angstrom
  a_info = {'colors':{'Ba':(.1,1,.1),'Cu':(.8,.8,0), 'As':(0,0,.7)}, 'radii':{'Ba':1.4}, 'bonds':{'Cu_As':5}}
  model = Model.Model(params=a_info, fname='data_files/BaCuAs.scf.in')
  xcp = XCP_Atoms.XCP_Atoms(model=model)
  xcp.start_crystal_view()

