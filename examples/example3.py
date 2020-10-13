from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # 2D SnTe with colors red & blue for Sn & Te respectively.
  # Bonds are drawn between Sn & Te atoms if they are closer together than 8 angstrom
  cpy = XCP.XtraCrysPy(inputfile='SnTe_data/SnTe.scf.in', species={'Sn':{'color':(1,0,0),'radius':.5}, 'Te':{'color':(0,0,1),'radius':.8}}, bonds={'Sn_Te':8.})

  cpy.start_cryspy_view(title='SnTe', nx=4, ny=4)
