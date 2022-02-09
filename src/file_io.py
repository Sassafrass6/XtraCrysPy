

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


def read_band_path_PAO ( fname ):
  '''
  '''
  import numpy as np

  tags = []
  npnts = []
  with open(fname, 'r') as f:
    ls = f.readline().split()
    while not len(ls) == 0:
      tags.append(ls[0])
      npnts.append(int(ls[1]))
      ls = f.readline().split()

  kcnt = 0
  ftags = []
  findex = [0]
  for i in range(len(tags)-1):
    if tags[i] == 'G':
      tags[i] = r'$\Gamma$'

    if npnts[i] == 0:
      ftags[-1] += '|' + tags[i]
    else:
      ftags.append(tags[i])
      findex.append(npnts[i] + findex[-1])
  ftags.append(tags[-1] if not tags[-1]=='G' else r'$\Gamma$')
    
  return findex, ftags


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


def read_transport_PAO ( fname ):

  import numpy as np

  nene = 0
  ntemp = 0
  enes = temps = tensors = None

  with open(fname, 'r') as f:
    lines = f.readlines()

    nl = len(lines)
    ftemp = float(lines[ntemp].split()[0])
    ptemp = ftemp
    while ftemp == ptemp and nene < nl-1:
      nene += 1
      ptemp = float(lines[nene].split()[0])
    if nene == nl-1:
      nene += 1

    while ntemp*nene < nl:
      ntemp += 1

    enes = np.empty(nene, dtype=float)
    temps = np.empty(ntemp, dtype=float)
    tensors = np.empty((ntemp,nene,3,3), dtype=float)

    iL = 0
    for i in range(ntemp):
      temps[i] = float(lines[iL].split()[0])
      for j in range(nene):
        ls = lines[iL].split()
        if i == 0:
          enes[j] = float(ls[1])
        for k in range(3):
          tensors[i,j,k,k] = float(ls[2+k])
        for ik,k in enumerate([(0,1), (0,2), [1,2]]):
          tensors[i,j,k[0],k[1]] = tensors[i,j,k[1],k[0]] = float(ls[5+ik])
        iL += 1
 
    return enes, temps, tensors


def read_transport_PAO_parse ( fname ):
  '''
  '''
  import numpy as np

  enes = []
  temps = [-1]
  tensors = []
  with open(fname, 'r') as f:
    tene = []
    ttemp = []
    ttensor = []
    for l in f.readlines():
      ls = l.split()
      ntemp = float(ls[0])
      if not temps[-1] == ntemp:
        enes.append(tene); tene = []
        temps.append(ntemp)
  temps.pop(0)



