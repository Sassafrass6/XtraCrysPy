

def read_dos_QE ( fname, return_ef=True ):
  '''
  '''
  import numpy as np

  ef = None
  es = []
  dos = []
  int_dos = []
  with open(fname, 'r') as f:
    lines = f.readlines()

    ef = float(lines[0].split()[8])

    for l in lines[1:]:
      ls = l.split()
      for i,a in enumerate([es,dos,int_dos]):
        a.append(float(ls[i]))

  es = np.array(es)
  dos = np.array(dos)
  int_dos = np.array(int_dos)

  return (ef, es, dos) if return_ef else (es, dos)


def read_bands_QE_dat ( fname ):
  '''
  '''
  import numpy as np

  ks = []
  bands = []
  with open(fname, 'r') as f:
    lines = f.readlines()

    ls = lines[0].split()
    nbnd = int(ls[2][:-1])
    nks = int(ls[4])

    li = 1
    while len(ks) < nks:
      ks.append([float(v) for v in lines[li].split()])
      band = []
      while len(band) < nbnd:
        li += 1
        band += [float(v) for v in lines[li].split()]
      bands.append(band)
      li += 1

  ks = np.array(ks)
  bands = np.array(bands).T

  return ks, bands


def read_bands_QE_agr ( fname ):
  '''
  '''
  import numpy as np

  bands = []
  with open(fname, 'r') as f:
    band = []
    for l in f.readlines():
      ls = l.split()
      if len(ls) == 0:
        bands.append(band)
        band = []
      else:
        band.append(float(ls[1]))
  return np.array(bands)


def read_dos_PAO ( fname ):
  '''
  '''
  import numpy as np

  es = []
  dos = []
  with open(fname, 'r') as f:
    for l in f.readlines():
      ls = l.split()
      es.append(float(ls[0]))
      dos.append(float(ls[1]))

  es = np.array(es)
  dos = np.array(dos)

  return es, dos


def read_bands_PAO ( fname ):
  '''
  '''
  import numpy as np

  bands = []
  with open(fname, 'r') as f:
    for l in f.readlines():
      bands.append([float(v) for v in l.split()[1:]])

  return np.array(bands).T
