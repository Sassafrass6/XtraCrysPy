from XCrysPy import XCrysPy as XCP

if '__main__' == __name__:

  # Define params
  a = 5.65

  # Unit origin
  origin = [0,0,0]

  # Cell parameters (Normalized)
  fcc_vecs = [[a/2,a/2,0], [0,a/2,a/2], [a/2,0,a/2]]

  # Atomic basis
  atoms = [[0,0,0], [a/4,a/4,a/4]]

  # Label for each atom and colors
  spec = ['Ga', 'As']
  spec_col = {'Ga':(1,0,0), 'As':(0,0,1)}

  # Plot with CrysPy
  cpy = XCP.XCrysPy(lattice=fcc_vecs, basis=atoms, species=spec, spec_col=spec_col, origin=origin, nx=3, ny=3, nz=3, bond_dists=3., boundary=False)

