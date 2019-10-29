from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron with modified background and boundary colors.
  #    Tuples formatted as (R,G,B)
  crystal = XCP.XtraCrysPy(inputfile='Fe.scf.in', bg_col=(1,1,1), bnd_col=(0,0,0), boundary=True)
