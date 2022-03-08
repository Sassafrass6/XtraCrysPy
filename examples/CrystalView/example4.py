from XtraCrysPy import XCP_Atoms as XCP

if '__main__' == __name__:

  # Relaxation
  a_info = {'colors':{'Ca':(0,1,0),'Al':(0,0,1),'P':(1,0,0)}, 'bonds':{'Al_P':5, 'Ca_Al':5, 'Ca_P':5}}
  xcp = XCP.XCP_Atoms(params=a_info, model='data_files/CaAlP.relax.out', relax=True, sel_type='Chain')
  xcp.start_crystal_view()

