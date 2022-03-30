from numpy import ndarray


def read_relaxed_coordinates_QE ( fname:str ):
  '''
    Reads relaxed atomic positions from a QE .out file. If vcrelax is set True, the crystal coordinates are also read.

    Arguments:
      fname (str): File name (including path) for the .out file
      vcrelax (bool): True reads crystal coordinates in addition to atomic positions
      read_all (bool): True forces all relax steps to be read. If EoF is encountered before 'final coordinates' the last coordinates to appear in the file are retunred. If no coordinates are found, an empty dictionary is returned.

    Returns:
      (dict): Dictionary with one or two entries - 'apos' for atomic positions and 'coord' for crystal coordinates.
  '''
  from os.path import isfile,join
  from os import getcwd
  import numpy as np

  if not isfile(fname):
    raise FileNotFoundError('File {} does not exist.'.format(join(getcwd(),fname)))

  abc = []
  ibrav = 0
  species = []
  cell_params = []
  celldm = np.empty(6, dtype=float)
  struct = {'lunit':'bohr', 'aunit':'alat'}
  with open(fname, 'r') as f:
    lines = f.readlines()

    eL = 0
    nL = len(lines)

    try:
      def read_apos ( sind ):
        apos = []
        while lines[sind] != '\n' and not 'End final coordinates' in lines[sind]:
          apos.append([float(v) for v in lines[sind].split()[1:4]])
          sind += 1
        return sind, apos

      while 'bravais-lattice' not in lines[eL]:
        eL += 1
      ibrav = int(lines[eL].split()[3])


      while  'celldm' not in lines[eL]:
        eL += 1
      celldm[:3] = [float(v) for i,v in enumerate(lines[eL].split()) if i%2==1]
      celldm[3:] = [float(v) for i,v in enumerate(lines[eL+1].split()) if i%2==1]

      if ibrav != 0:
        from .lattice_format import lattice_format_QE
        cell_params.append(lattice_format_QE(ibrav, celldm))
      else:
        while 'crystal axes' not in lines[eL]:
          eL += 1
        coord = []
        for l in lines[eL+1:eL+4]:
          coord.append([celldm[0]*float(v) for v in l.split()[3:6]])
        cell_params.append(coord)

      while 'site n.' not in lines[eL]:
        eL += 1
      eL += 1
      apos = []
      while 'End' not in lines[eL] and lines[eL] != '\n':
        line = lines[eL].split()
        species.append(line[1])
        apos.append([float(v) for v in line[6:9]])
        eL += 1
      abc.append(apos)

      while eL < nL:
        while eL < nL and 'CELL_PARAMETERS' not in lines[eL] and 'ATOMIC_POSITIONS' not in lines[eL]:
          eL += 1
        if eL >= nL:
          break
        if 'ATOMIC_POSITIONS' in lines[eL]:
          unit = lines[eL].split()[1].strip('(){{}}')
          if len(unit) > 1:
            struct['aunit'] = unit
          eL,apos = read_apos(eL+1)
          abc.append(apos)
        elif 'CELL_PARAMETERS' in lines[eL]:
          coord = []
          unit = lines[eL].split()[1].strip('(){{}}')
          if len(unit) > 1:
            struct['lunit'] = unit
          for l in lines[eL+1:eL+4]:
            coord.append(np.array([float(v) for v in l.split()]))
          cell_params.append(coord)
          eL += 4

          while 'ATOMIC_POSITIONS' not in lines[eL]:
            eL += 1
          eL,apos = read_apos(eL+1)
          abc.append(apos)

    except Exception as e:
      print('WARNING: No atomic positions or cell coordinates were found.', flush=True)
      raise e

  struct['species'] = species
  struct['lattice'] = np.array(cell_params)
  struct['abc'] = np.array(abc)

  return struct


def struct_from_inputfile_QE ( fname:str ) -> dict:
  '''
    Generate a dictionary containing all atomic information from a QE inputfile
    WARNING: Currently only the control blocks are read. Atomic cards are not...

    Arguments:
      fname (str): Name (including path) of the inputfile

    Returns:
      (dict): Structure dictionary
  '''
  from os.path import isfile
  import numpy as np
  import re

  if not isfile(fname):
    raise FileNotFoundError('File {} does not exist.'.format(fname))

  pattern = re.compile('&(.*?)/@')
  nocomma = lambda s : s.replace(',', '')
  inquote = lambda s : s.split('\'')[1]
  qenum = lambda s : s.split('=')[1].replace('d', 'e')
  qeint = lambda s : int(qenum(s))
  qefloat = lambda s : float(qenum(s))
  qebool = lambda s : True if s.split('.')[1][0].lower() == 't' else False

  struct = {}
  celldm = np.zeros(6, dtype=float)
  with open(fname, 'r') as f:
    fstr = f.read()
    matches = pattern.findall(fstr.replace(' ', '').replace('\n', '@ '))
    for m in matches:
      m = [s.replace(' ', '') for s in re.split(', |@', m)]
      mf = []
      for v in m:
        mf += v.split(',')
      block = mf.pop(0).lower()
      if 'control' in block:
        for s in mf:
          if 'calculation' in s:
            struct['calc'] = inquote(s)
          elif 'outdir' in s:
            struct['outdir'] = inquote(s)
          elif 'prefix' in s:
            struct['id'] = inquote(s)
          elif 'pseudo_dir' in s:
            struct['ppdir'] = inquote(s)
          elif 'tstress' in s:
            struct['tstress'] = qebool(s)
          elif 'tprnfor' in s:
            struct['tprnfor'] = qebool(s)
      elif 'system' in block:
        for s in mf:
          if 'ibrav' in s:
            struct['ibrav'] = qeint(s)
          if 'nat' in s:
            struct['nat'] = qeint(s)
          elif 'ecutwfc' in s:
            struct['ecutwfc'] = qefloat(s)
          elif 'ecutrho' in s:
            struct['ecutrho'] = qefloat(s)
          elif 'occupations' in s:
            struct['occ'] = inquote(s)
          elif 'smearing' in s:
            struct['smearing'] = inquote(s)
          elif 'degauss' in s:
            struct['degauss'] = qefloat(s)
          elif 'celldm' in s:
            cpattern = re.compile('\(([^\)]+)\)')
            dm = int(cpattern.findall(s)[0]) - 1
            celldm[dm] = qefloat(s)
      elif 'electrons' in block:
        for s in mf:
          if 'conv_thr' in s:
            struct['conv_thr'] = qefloat(s)
      elif 'ions' in block:
        pass
      elif 'cell' in block:
        pass

  if struct['ibrav'] != 0:
    struct['lunit'] = 'bohr'
    struct['celldm'] = celldm

  il = 0
  cl = 0
  kl = 0
  fstr = fstr.split('\n')
  while not 'ATOMIC_POSITION' in fstr[il]:
    if 'CELL_PARAM' in fstr[il]:
      cl = il
    if 'K_POINTS' in fstr[il]:
      kl = il
    il += 1
  unit = fstr[il].split()
  struct['aunit'] = unit[1].strip('(){{}}') if len(unit) > 1 else 'alat'
  il += 1

  spec = []
  abc = np.empty((struct['nat'],3), dtype=float)
  for i in range(struct['nat']):
    ls = fstr[il+i].split()
    spec.append(ls[0])
    abc[i,:] = np.array([float(v) for v in ls[1:]])

  struct['abc'] = abc
  struct['species'] = spec

  if kl == 0:
    while not 'K_POINTS' in fstr[kl]:
      if 'CELL_PARAM' in fstr[cl]:
        cl = kl
      kl += 1

  if 'gamma' in fstr[kl]:
    struct['kpnts'] = {'type':'gamma'}
  elif 'automatic' in fstr[kl]:
    nk1,nk2,nk3,o1 = (int(v) for v in fstr[kl+1].split()[:4])
    struct['kpnts'] = {'type':'automatic', 'nk1':nk1, 'nk2':nk2, 'nk3':nk3, 'offset':o1}

  nf = len(fstr)
  if cl == 0:
    while cl < nf and not 'CELL_PARAM' in fstr[cl]:
      cl += 1

  if cl < nf:
    unit = fstr[cl].split()
    struct['lunit'] = unit[1].strip('(){{}}') if len(unit) > 1 else 'alat'
    cl += 1
    struct['lattice'] = np.array([np.array([float(v) for v in fstr[cl+c].split()]) for c in range(3)])

  if 'lattice' not in struct:
    from .lattice_format import lattice_format_QE
    struct['lattice'] = lattice_format_QE(struct['ibrav'], struct['celldm'])

  if struct['lunit'] == 'angstrom':
    struct['lattice'] *= 1.88973
  elif struct['lunit'] != 'bohr':
    print('WARNING: Unit {} may behave strangely'.format(struct['lunit']))

  if struct['aunit'] == 'bohr' or struct['aunit'] == 'angstrom':
    if struct['aunit'] == 'angstrom':
      struct['abc'] *= 1.88973
    struct['abc'] = struct['abc'] @ np.linalg.inv(struct['lattice'])

  return struct


def struct_from_inputfile_CIF ( fname:str ) -> dict:
  '''
  '''
  import numpy as np
  import re

  def strip_float ( text:str ):
    pat = '[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)'
    return float(re.search(pat, text).group())

  lines = []
  struct = {}
  with open(fname, 'r') as f:
    lines = ' '.join([s.replace('\n', ' ') for s in f.readlines() if not s[0]=='#']).split()

  i = 0
  flines = []
  while i < len(lines):
    c = None
    l = lines[i]
    if l[0] == "'" or l[0] == '"' or l[0] == ';':
      c = l[0]

    i += 1
    if c is None:
      flines.append(lines[i-1])

    else:
      fstr = [l]
      while lines[i][-1] != c:
        fstr.append(lines[i])
        i += 1
      fstr.append(lines[i])
      flines.append(' '.join(fstr))
      i += 1
  lines = flines

  eL = 0
  nL = len(lines)
  data_blocks = {}
  while eL < nL:
    while eL < nL and 'data_' not in lines[eL]:
      eL += 1
    if eL == nL:
      break

    
    eL += 1
    data = {}
    loops = []
    prefix = '_'.join(lines[eL-1].split('_')[1:])
    data_blocks[prefix] = {'data':data, 'loops':loops}
    while eL < nL and (lines[eL][0] == '_' or 'loop_' in lines[eL] ):

      tag = lines[eL]
      if lines[eL][0] == '_':
        data[tag] = lines[eL+1]
        eL += 2
      elif 'loop_' in tag:
        nele = 0
        labels = []
        while lines[eL+nele+1][0] == '_':
          nele += 1
          labels.append(lines[eL+nele])
        eles = []
        eL += nele + 1
        while eL < nL and lines[eL][0] != '_' and 'loop_' not in lines[eL]:
          eles.append(lines[eL:eL+nele])
          eL += nele
        loops.append((labels,eles))
      else:
        raise Exception('Unknown field {}'.format(tag))
      '''
      continue
        #  ind = ord(tag.split('_')[-1]) - 97
        #  #cell_param[ind] = strip_float(lines[eL+1])
        #  data[tag] = strip_float(lines[eL+1])
      tag = lambda t : t.split('_')[-1]
      elif '_cell_angle_' in lines[eL]:
        ind = 3 + ['alpha', 'beta', 'gamma'].index(tag(ls[0]))
      '''

  print(data_blocks)
  quit()
    

def struct_from_inputfile_POSCAR ( fname:str ) -> dict:
  import numpy as np

  struct = {}
  with open(fname, 'r') as f:

    lines = f.readlines()

    alat = float(lines[1].split()[0])
    lattice = np.empty((3,3), dtype=float)
    for i in range(3):
      lattice[i,:] = np.array([float(v) for v in lines[2+i].split()])

    uspec = lines[5].split()
    nspec = lines[6].split()
    numeric = nspec[0].isnumeric()
    ucount = [int(v) for v in (nspec if numeric else uspec)]

    spec = []
    for i,c in enumerate(ucount):
      spec += c * [(uspec[i] if numeric else str(i))]

    lc = 7 if numeric else 6
    if lines[lc][0].lower() == 's':
      lc += 1

    cartesian = lines[lc][0].lower() in ['c', 'k']

    lc += 1
    nat = len(spec)
    abc = np.empty((nat,3), dtype=float)
    for i in range(nat):
      abc[i,:] = np.array([float(v) for v in lines[lc+i].split()])

    lattice *= 1.88973 * alat
    if cartesian:
      abc = 1.88973 * abc @ np.linalg.inv(lattice)

    struct['abc'] = abc
    struct['species'] = spec
    struct['lattice'] = lattice

    return struct


def infer_file_type ( fname:str ):

  extension = fname.split('.')
  if len(extension) <= 1:
    print('Cannot infer file type without an extension. Assuming qe')
    return 'qe'

  extension = extension[-1].lower()
  if extension in ['in', 'out']:
    ftype = 'qe'
  elif extension == 'poscar':
    ftype = 'poscar'
  else:
    ftype = extension

  return ftype


def read_relaxed_coordinates ( fname:str, ftype='automatic' ):
  '''
    Assumed that the file type is a QE relax output file
  '''

  if ftype == 'automatic':
    ftype = infer_file_type(fname)

  if ftype == 'qe':
    return read_relaxed_coordinates_QE(fname)
  else:
    raise ValueError('Cannot read relax file type {}'.format(ftype))


def struct_from_inputfile ( fname:str, ftype='automatic' ):
  '''
  '''

  if ftype == 'automatic':
    ftype = infer_file_type(fname)

  try:
    if ftype == 'qe':
      return struct_from_inputfile_QE(fname)
    elif ftype == 'cif':
      return struct_from_inputfile_CIF(fname)
    elif ftype == 'poscar':
      return struct_from_inputfile_POSCAR(fname)
    else:
      raise ValueError('Cannot read file type {}'.format(ftype))

  except Exception as e:
    print('Failed to read file {}'.format(fname))
    raise e


def read_dos_QE ( fname:str, read_ef=True ):
  '''
  '''
  import numpy as np

  ef = None
  es = []
  dos = []
  int_dos = []
  with open(fname, 'r') as f:
    lines = f.readlines()

    ef = float(lines[0].split()[8]) if read_ef else 0.

    for l in lines[1:]:
      ls = l.split()
      for i,a in enumerate([es,dos,int_dos]):
        a.append(float(ls[i]))

  es = np.array(es)
  dos = np.array(dos)
  int_dos = np.array(int_dos)

  return (ef, es, dos)


def read_bands_QE_dat ( fname:str ):
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


def read_bands_QE_agr ( fname:str ):
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


def read_band_path_PAO ( fname:str ):
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
    if tags[i] == 'G' or tags[i] == 'gG':
      tags[i] = r'$\Gamma$'

    if npnts[i] == 0:
      ftags[-1] += '|' + tags[i]
    else:
      ftags.append(tags[i])
      findex.append(npnts[i] + findex[-1])
  ftags.append(tags[-1] if not tags[-1]=='G' else r'$\Gamma$')
    
  return findex, ftags


def read_dos_PAO ( fname:str ):
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


def read_bands_PAO ( fname:str ):
  '''
  '''
  import numpy as np

  bands = []
  with open(fname, 'r') as f:
    for l in f.readlines():
      bands.append([float(v) for v in l.split()[1:]])

  return np.array(bands).T


def read_transport_PAO ( fname:str ):

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


def read_datablocks_XSF ( fname:str ):
  '''
  '''
  import numpy as np

  blocks = []
  with open(fname, 'r') as f:
    lines = f.readlines()
    nlines = len(lines)

    ind = 0
    while ind < nlines:
      while 'BEGIN_DATAGRID' not in lines[ind]:
        ind += 1
        if ind == nlines:
          return np.array(blocks)
      if 'BEGIN_DATAGRID' in lines[ind]:
        ind += 1
        dims = np.array([int(v) for v in lines[ind].split()])
        ind += 5
        data = np.empty(dims, dtype=float)
        for i in range(dims[0]):
          for j in range(dims[1]):
            data[i,j,:] = np.array([float(v) for v in lines[ind].split()])
            ind += 1
        blocks.append(data)
      ind += 1

  return np.array(blocks)
