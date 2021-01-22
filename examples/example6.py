from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron with modified background and boundary colors.
  #    Tuples formatted as (R,G,B)
  species = {'Si':{'color':(0,.8,.1)}}
  bonds = {'Si_Si':2.4*1.88973}

  crystal = XCP.XtraCrysPy(inputfile='Si.vcrelax.out', relax=True, species=species, bonds=bonds)
  crystal.start_cryspy_view(title='Si Relax', bg_color=(0,0,.1), f_color=(.6,.2,0), nx=2, nz=2)
