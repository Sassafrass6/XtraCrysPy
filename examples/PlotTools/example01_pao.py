
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot,QE_Plot

  f_dos = './output_PAO/dosdk_0.dat'
  f_band = './output_PAO/bands_0.dat'
  f_symp = './output_PAO/kpath_points.txt'
  f_sigma = './output_PAO/sigmadk_0.dat'
  f_seebeck = './output_PAO/Seebeck_0.dat'

  pplt = PAO_Plot.PAO_Plot()

  pplt.plot_dos_PAO(f_dos, title='Si2 FCC DoS', vertical=False)
  pplt.plot_bands_PAO(f_band, sym_points=f_symp)
  pplt.plot_dos_beside_bands_PAO(f_dos, f_band, sym_points=f_symp, y_lim=(-11.5,5), dos_ticks=True)

  pplt.plot_electrical_conductivity_PAO(f_sigma, title='Ba8Cu16As30 Conductivity', x_lim=(-5,1))
  pplt.plot_seebeck_PAO(f_seebeck, t_ele=[], x_lim=(-.4,.8), y_lim=(-1.1e3,1.1e3), col='black')

