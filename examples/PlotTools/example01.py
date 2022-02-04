
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot,QE_Plot

  qplt = QE_Plot.QE_Plot()
  qplt.plot_bands_QE('./output_QE/BaCuAs.agr')
  qplt.plot_dos_beside_bands_QE('./output_QE/BaCuAs.dos', './output_QE/BaCuAs.agr')
  quit()

  pplt = PAO_Plot.PAO_Plot()
  pplt.plot_bands_PAO('./output_PAO/bands_0.dat', y_lim=(-11.5,5))
  pplt.plot_dos_beside_bands_PAO('./output_PAO/dosdk_0.dat', './output_PAO/bands_0.dat', y_lim=(-11.5,5), dos_ticks=True)

  sym_point_positions = [0,7,17,24,34,46,53,65,72,79,86]
  sym_point_tags = ['$\Gamma$','Z','D','B','$\Gamma$','A','E','Z','C2','Y2','$\Gamma$']
  sym_points = (sym_point_positions, sym_point_tags)
  pplt.plot_dos_beside_bands_PAO('./output_PAO/dosdk_0.dat', './output_PAO/bands_0.dat', sym_points=sym_points)

