

class QE_Plot:

  def __init__ ( self ):
    pass

  def plot_dos_QE ( self, fname, title=None, x_lim=None, y_lim=None, vertical=False, col='black' ):
    '''
    '''
    from .file_io import read_dos_QE
    from .plot_functions import plot_dos

    ef, es, dos = read_dos_QE(fname)
    plot_dos(es-ef, dos, title, x_lim, y_lim, vertical, col)


  def plot_bands_QE ( self, fname, eF=0., sym_points=None, title=None, y_lim=None, col='black' ):
    '''
    '''
    from .plot_functions import plot_bands

    bands = None
    ext = fname.split('.')[-1]
    if ext == 'agr':
      from .file_io import read_bands_QE_agr
      bands = read_bands_QE_agr(fname)
    elif ext == 'dat':
      from .file_io import read_bands_QE_dat
      ks,bands = read_bands_QE_dat(fname)
    else:
      raise Exception('Cannot read band file with extension .{}'.format(ext))

    if type(eF) is str:
      from .file_io import read_dos_QE
      eF,_,_ = read_dos_QE(eF)

    plot_bands(bands-eF, sym_points, title, None, y_lim, col)


  def plot_dos_beside_bands_QE ( self, fn_dos, fn_bands, align_eF=True, sym_points=None, title=None, x_lim=None, y_lim=None, col='black', dos_ticks=False ):
    '''
    '''
    from .file_io import read_dos_QE
    from .plot_functions import plot_dos_beside_bands

    bands = None
    ext = fn_bands.split('.')[-1]
    if ext == 'agr':
      from .file_io import read_bands_QE_agr
      bands = read_bands_QE_agr(fn_bands)
    elif ext == 'dat':
      from .file_io import read_bands_QE_dat
      ks, bands = read_bands_QE_dat(fn_bands)
    else:
      raise Exception('Cannot read band file with extension .{}'.format(ext))

    ef, es, dos = read_dos_QE(fn_dos)

    if align_eF:
      es -= ef
      bands -= ef

    plot_dos_beside_bands(es, dos, bands, sym_points, title, None, x_lim, y_lim, col, dos_ticks)
