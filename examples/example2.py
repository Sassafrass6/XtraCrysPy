from XCrysPy import XCrysPy as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron
  crystal = XCP.XCrysPy(qe_fname='Fe.scf.in')
  crystal.draw_cell(boundary=True)
