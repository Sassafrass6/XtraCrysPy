from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron with modified background and boundary colors.
  #    Tuples formatted as (R,G,B)
  species = {'Ca':{'color':(0,1,0)}, 'Al':{'color':(0,0,1)}, 'P':{'color':(1,0,0)}}
  bonds = {'Al_P':2.6*1.88973}

  crystal = XCP.XtraCrysPy(inputfile='CaAlP.relax.out', relax=True, species=species, bonds=bonds)
  crystal.start_cryspy_view(title='Clathrate Relax', bg_color=(0,0,.1), f_color=(.6,.2,0), nx=1, nz=1)
