from XtraCrysPy import XCP_Atoms as XCP

if '__main__' == __name__:

  # Draw the unit cell for Iron with modified background and boundary colors.
  #    Tuples formatted as (R,G,B)
#  species = {'Ca':{'color':(0,1,0)}, 'Al':{'color':(0,0,1)}, 'P':{'color':(1,0,0)}}
#  bonds = {'Al_P':2.6*1.88973}

  a_info = {'colors':{'Ca':(0,1,0),'Al':(0,0,1),'P':(1,0,0)}, 'bonds':{'Al_P':5, 'Ca_Al':5, 'Ca_P':5}}
  xcp = XCP.XCP_Atoms(params=a_info, model='CaAlP.relax.out', relax=True, sel_type='chain')
  xcp.start_crystal_view()

