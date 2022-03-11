from XtraCrysPy import XCP_Atoms as XCP

if '__main__' == __name__:

  a_info = {'bonds':{'In_As':6.1, 'As_As':6.1, 'In_In':6.1}}
  xcp = XCP.XCP_Atoms(params=a_info, model='data_files/CsInAs.scf.in', sel_type='Chain')
  xcp.start_crystal_view()

