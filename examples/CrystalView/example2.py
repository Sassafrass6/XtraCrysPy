from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron with modified background and boundary colors.
  #    Tuples formatted as (R,G,B)
  crystal = XCP.XtraCrysPy(inputfile='Fe.scf.in')
  crystal.start_cryspy_view(title='Fe', bg_color=(1,1,1), f_color=(0,0,0))
