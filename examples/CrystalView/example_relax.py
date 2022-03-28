from XtraCrysPy import XCP_Atoms as XCP

if '__main__' == __name__:

  a_info = {'colors':{'Ca':(0,1,0),'Al':(0,0,1),'P':(1,0,0)}, 'bonds':5}
  # Relaxation
  xcp = XCP.XCP_Atoms(model='data_files/CaAlP.relax.out', 
                      params=a_info, relax=True, sel_type='Chain')
  xcp.start_crystal_view()

