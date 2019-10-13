from XCrysPy import XCrysPy

if '__main__' == __name__:

  # Draw the unit cell for Iron
  crystal = XCrysPy(qe_fname='Fe.scf.in')
  crystal.draw_cell(boundary=True)
