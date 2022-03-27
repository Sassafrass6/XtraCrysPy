
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot

  f_symp = './output_PAO/Fe.kpath_points.txt'
  f_band = './output_PAO/Fe.bands_0.dat'
  f_berry = './output_PAO/Fe.Omega_z_xy.dat'

  pplt = PAO_Plot.PAO_Plot()

  # Functions arguments (tiltle y_lim, etc) can be used in any of the plot functions
  pplt.plot_berry(f_berry, sym_points=f_symp)
  pplt.plot_bands(f_band, sym_points=f_symp, y_lim=(-4,4))
  berry_label = '$\Omega^{z}$($\\bf{k}$)'
  pplt.plot_berry_under_bands(f_berry, f_band, sym_points=f_symp, y_lim=(-4,4), dos_ticks=True, berry_label=berry_label)
