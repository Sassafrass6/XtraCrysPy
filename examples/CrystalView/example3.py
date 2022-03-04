from XtraCrysPy import XtraCrysPy,Model


if '__main__' == __name__:

  # 2D SnTe with colors red & blue for Sn & Te respectively.
  # Bonds are drawn between Sn & Te atoms if they are closer together than 5 angstrom
  a_info = {'colors':{'Sn':(1,0,0),'Te':(0,0,1)}, 'radii':{'Sn':1,'Te':1.1}, 'bonds':{'Sn_Te':6}}
  model = Model.Model(params=a_info, fname='SnTe_data/SnTe.scf.in')
  xcp = XtraCrysPy.XtraCrysPy()
  xcp.render_atomic_model(model, nsc=(3,3,1))#, bond_type='primary')
  xcp.start_crystal_view()

