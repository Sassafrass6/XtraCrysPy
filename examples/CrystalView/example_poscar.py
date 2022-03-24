from XtraCrysPy import XCP_Atoms,Model

if '__main__' == __name__:

  fname = 'data_files/BN.poscar'
  xcp = XCP_Atoms.XCP_Atoms(model=fname, nsc=(2,2,2), perspective=True)
  xcp.start_crystal_view()

