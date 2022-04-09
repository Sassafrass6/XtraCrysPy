from XtraCrysPy.Atomic import Atomic

if '__main__' == __name__:

  a_info = {'bonds':4.2}
  fname = 'data_files/BaO3Ti.cif'
  xcp = Atomic(params=a_info, model=fname, sel_type='Distance', nsc=(2,2,2))
  xcp.start_crystal_view()

