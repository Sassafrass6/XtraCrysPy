from matplotlib import pyplot as plt

def plot_dos ( fname, title=None, x_lim=None, y_lim=None, vertical=False, col='black' ):
  tit = 'DoS' if title is None else title
  fig.suptitle(tit)
  pass

def plot_bands ( bands, sym_points, title, y_lim, col ):
  '''
  '''

  fig = plt.figure()
  ax = fig.add_subplot(111)

  tit = 'Band Structure' if title is None else title
  fig.suptitle(tit)

  for b in bands:
    ax.plot(b, color=col)
  if y_lim is None:
    y_lim = ax.get_ylim()
  ax.set_ylim(*y_lim)
  if sym_points is None:
    ax.xaxis.set_visible(False)
  else:
    ax.set_xticks(sym_points[0])
    ax.set_xticklabels(sym_points[1])
    ax.vlines(sym_points[0], y_lim[0], y_lim[1], color='gray')
  ax.set_ylabel('Energy (eV)', fontsize=12)

  plt.show()


def plot_dos_beside_bands ( es, dos, bands, sym_points, title, x_lim, y_lim, col, dos_ticks ):
  '''
  '''
  from matplotlib import gridspec

  fig = plt.figure()
  spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[5,1])

  tit = 'Band Structure and DoS' if title is None else title
  fig.suptitle(tit)

  ax_b = fig.add_subplot(spec[0])
  ax_d = fig.add_subplot(spec[1])

  for b in bands:
    ax_b.plot(b, color=col)
  if y_lim is None:
    y_lim = ax_b.get_ylim()
  ax_b.set_ylim(*y_lim)
  if sym_points is None:
    ax_b.xaxis.set_visible(False)
  else:
    ax_b.set_xticks(sym_points[0])
    ax_b.set_xticklabels(sym_points[1])
    ax_b.vlines(sym_points[0], y_lim[0], y_lim[1], color='gray')
  ax_b.set_ylabel('Energy (eV)', fontsize=12)
  
  ax_d.plot(dos, es, color=col)
  if not x_lim is None:
    ax_d.set_xlim(*x_lim)
  if not y_lim is None:
    ax_d.set_ylim(*y_lim)
  if not dos_ticks:
    ax_d.yaxis.set_visible(False)
    plt.tight_layout()

  plt.show()


def plot_electrical_conductivity ( fname, t_ele, x_lim, y_lim, col, vT ):
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

def plot_seebeck ( x_lim, y_lim, col ):
  pass
