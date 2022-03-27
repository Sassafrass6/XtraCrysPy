
if __name__ == '__main__':
  from XtraCrysPy import PAO_Plot
  from glob import glob

  pplt = PAO_Plot.PAO_Plot()

  # Plot a single shc file
  f_shc = './output_PAO/Pt.shcEf_z_xy.dat'
  pplt.plot_shc(f_shc, title='Pt SHC', legend=False, cols='black')

  # Plot a list of shc files
  #   cols argument can also be provided as a list of strings or 3-colors
  files = glob('./output_PAO/Pt.shc*.dat')
  pplt.plot_shc(files, title='Pt SHC')

