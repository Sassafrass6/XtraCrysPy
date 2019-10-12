from CrysPy import CrysPy

if '__main__' == __name__:

  pyc = CrysPy(qe_fname='SnTe.scf.in', spec_col={'Sn':(1,0,0), 'Te':(0,0,1)})
  pyc.draw_cell(nx=3, ny=3, nz=1, boundary=False)
  pyc.draw_bonds(dist={'Sn_Te':8.})
