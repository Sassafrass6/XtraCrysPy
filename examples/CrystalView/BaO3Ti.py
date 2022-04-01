from XtraCrysPy import XCP_Atoms

if '__main__' == __name__:

  a_info = {'bonds':4.5}
  xcp = XCP_Atoms.XCP_Atoms(model='data_files/BaO3Ti.cif', params=a_info)
  xcp.start_crystal_view()

