from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  cpy = XCP.XtraCrysPy(perspective=False)
  cpy.plot_bxsf(fname='Spin_Berry_z_xy.bxsf', iso=[250])
