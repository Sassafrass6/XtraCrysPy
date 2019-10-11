from CrysPy import CrysPy

if '__main__' == __name__:

  # Define params
  a = 5.

  # Unit origin
  origin = [0,0,0]

  # Cell parameters (Normalized)
  fcc_vecs = [[a/2,a/2,0], [0,a/2,a/2], [a/2,0,a/2]]

  # Atomic basis
  atoms = [[0,0,0], [.5,.5,.5]]

  # Label for each atom and colors
  spec = ['Sn', 'Te']
  spec_col = {'Sn':(1,0,0), 'Te':(0,0,1)}


  pyc = CrysPy(lattice=fcc_vecs, basis=atoms, species=spec, spec_col=spec_col, origin=origin)
  pyc.draw_cell(nx=4, ny=4, nz=4)
  pyc.draw_bonds(dist=3)
