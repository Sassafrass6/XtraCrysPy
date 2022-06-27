from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  a_info = {'bonds':3}
  fname = 'data_files/Si46.scf.in'
  xcp = Atomic(params=a_info, model=fname,
               constrain_atoms=True, sel_type='Chain')
  xcp.start_crystal_view()

