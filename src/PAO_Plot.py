from matplotlib import pyplot as plt

class PAO_Plot:

  def __init__ ( self ):
    print('Should this really be a class?')

  def plot_dos_PAO ( self, fname, title=None, x_lim=None, y_lim=None, vertical=False, col='black' ):
    #tit = 'DoS' if title is None else title
    #fig.suptitle(tit)
    pass

  def plot_bands_PAO ( self, fname, sym_points=None, title=None, y_lim=None, col='black' ):
    '''
    '''
    from .file_io import read_bands_PAO
    from .plot_functions import plot_bands

    plot_bands(read_bands_PAO(fname), sym_points, title, y_lim, col)


  def plot_dos_beside_bands_PAO ( self, fn_dos, fn_bands, sym_points=None, title=None, x_lim=None, y_lim=None, col='black', dos_ticks=False ):
    '''
    '''
    from .file_io import read_dos_PAO,read_bands_PAO
    from .plot_functions import plot_dos_beside_bands

    bands = read_bands_PAO(fn_bands)
    es, dos = read_dos_PAO(fn_dos)

    plot_dos_beside_bands(es, dos, bands, sym_points, title, x_lim, y_lim, col, dos_ticks)    


  def plot_electrical_conductivity_PAO ( self, fname, t_ele=[(0,0),(1,1),(2,2)], x_lim=None, y_lim=None, col=[(1,0,0),(0,1,0),(0,0,1)], vT=None ):
    '''
      Plot the electrical conductivity.

      Arguments:
        fname (str): File name (including relative path)
        t_ele (list): Tensor elements as tuple pairs (e.g. (1,2) for (y,z)). Default behavior is to plot the 3 diagonal elements seprately. Providing an empty list will average the diagonal components
        x_lim (tuple): Pair of axis limits (x_min, x_max)
        y_lim (tuple): Pair of axis limits (y_min, y_max)
        col (list): A list of 3-tuples (R,G,B), one for each tensor element.
        vT (float): Set to an energy to plot the conductivity vs temperature. The value of conductivity is taken at the provided energy for each temperature.
    '''
    pass

  def plot_seebeck_PAO ( self, fname, x_lim=None, y_lim=None, col='black' ):
    pass
