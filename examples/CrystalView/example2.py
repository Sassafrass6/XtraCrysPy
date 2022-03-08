from XtraCrysPy import XCP_BZ,Model

if '__main__' == __name__:

  # BZ for FCC Iron, from a QE inputfile
  model = Model.Model(fname='data_files/Fe.scf.in')
  xcp = XCP_BZ.XCP_BZ(model=model)
  xcp.start_crystal_view()

