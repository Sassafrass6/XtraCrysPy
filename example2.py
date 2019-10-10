from CrysPy import CrysPy

if '__main__' == __name__:

  specD = {'Sn':(1,0,0), 'Te':(0,0,1)}
  pyc = CrysPy(qe_fname='SnTe.scf.in', spec_col=specD)
  pyc.draw_cell(nx=3, ny=3, nz=1, boundary=False)
  pyc.draw_bonds(dist={'Sn_Te':8.})

