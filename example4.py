from CrysPy import CrysPy

if '__main__' == __name__:

  # Define params
  a = 5.65
  fcc_vecs = [[a/2.,a/2.,0], [0,a/2.,a/2.], [a/2.,0,a/2.]]
  atoms = [[0,0,0]]

  pyc = CrysPy(lattice=fcc_vecs, basis=atoms)
  pyc.draw_brillouin_zone()
