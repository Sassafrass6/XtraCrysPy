from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  fname = 'data_files/BN.poscar'
  xcp = Atomic(model=fname, params={'bonds':1.6}, nsc=(2,2,2), perspective=True)
  xcp.start_crystal_view()

