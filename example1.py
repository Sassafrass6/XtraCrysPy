from CrysPy import CrysPy

if '__main__' == __name__:

  # Define params
  a = 5
  b = 4
  c = 5

  # Unit origin
  origin = [0,0,0]

  # Cell parameters (Normalized)
  fcc_vecs = [[.5,.5,0], [0,.5,.5], [.5,0,.5]]

  # Atomic basis, labels and colors
  atoms = [[0,0,0], [.5,.5,.5]]
  spec = ['Sn', 'Te']
  spec_col = {'Sn':(1,0,0), 'Te':(0,0,1)}


  pyc = CrysPy(lattice=fcc_vecs, basis=atoms, A=a, B=b, C=c, species=spec, spec_col=spec_col, orig=origin)
  pyc.draw_cell(nx=4, ny=4, nz=4, boundary=True)
  pyc.draw_bonds(dist=3)

