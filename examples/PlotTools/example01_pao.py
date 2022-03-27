
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot,QE_Plot

  f_dos = './output_PAO/Si.dosdk_0.dat'
  f_band = './output_PAO/Si.bands_0.dat'
  f_symp = './output_PAO/Si.kpath_points.txt'
  f_sigma = './output_PAO/Al.sigmadk_0.dat'
  f_seebeck = './output_PAO/AlP.Seebeck_0.dat'

  pplt = PAO_Plot.PAO_Plot()

  # Functions arguments (tiltle y_lim, etc) can be used in any of the plot functions
  pplt.plot_dos(f_dos, title='Si2 FCC DoS', vertical=False)
  pplt.plot_bands(f_band, sym_points=f_symp)
  pplt.plot_dos_beside_bands(f_dos, f_band, sym_points=f_symp, y_lim=(-11.5,5), dos_ticks=True)

  # Argument t_ele default is [[0,0], [1,1], [2,2]]
  #   which plots the 3 diagonal elements
  pplt.plot_electrical_conductivity(f_sigma, title='Al Conductivity', x_lim=(-5,1))

  # Setting t_ele to [] causes the diagonal elements to be averaged.
  pplt.plot_seebeck(f_seebeck, title='AlP Seebeck', t_ele=[], col='black', x_lim=(-1,4), y_lim=(-2e4,2e4))

