from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  fname = 'data_files/CsInAs.cp2k'
  xcp = Atomic(model=fname, params={'bonds':3.3})
  xcp.start_crystal_view()

