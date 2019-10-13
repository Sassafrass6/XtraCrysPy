from CrysPy import CrysPy

if '__main__' == __name__:

  # Draw the unit cell for Iron
  crystal = CrysPy(qe_fname='Fe.scf.in')
  crystal.draw_cell(boundary=True)
