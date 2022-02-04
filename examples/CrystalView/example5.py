from XtraCrysPy import XtraCrysPy as XCP

if '__main__' == __name__:

  cpy = XCP.XtraCrysPy()
  cpy.plot_bxsf(fname='Spin_Berry_z_xy.bxsf', iso=[150], bands=[1], colors=[(0,1,0)], normals=True)
