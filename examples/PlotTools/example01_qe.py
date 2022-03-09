
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot,QE_Plot

  f_dos = './output_QE/BaCuAs.dos'
  f_band = './output_QE/BaCuAs.agr'

  # Functions arguments (y_lim, sym_points, etc) can be used in any of the plot functions
  qplt = QE_Plot.QE_Plot()
  qplt.plot_dos_QE(f_dos, vertical=True)

  sym_point_positions = [0,7,17,24,34,46,53,65,72,79,86]
  sym_point_tags = ['$\Gamma$','Z','D','B','$\Gamma$','A','E','Z','C2','Y2','$\Gamma$']
  sym_arg = (sym_point_positions, sym_point_tags)

  qplt.plot_bands_QE(f_band, y_lim=(8.5,10))
  qplt.plot_dos_beside_bands_QE(f_dos, f_band, sym_points=sym_arg)
