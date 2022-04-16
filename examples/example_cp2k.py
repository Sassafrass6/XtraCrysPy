from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  fname = 'data_files/CsInAs.cp2k'
  xcp = Atomic(model=fname, params={'bonds':6.1})
  xcp.start_crystal_view()

