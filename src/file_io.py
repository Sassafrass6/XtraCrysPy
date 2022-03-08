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
          eL,apos = read_apos(eL+1)
          abc.append(apos)
        elif 'CELL_PARAMETERS' in lines[eL]:
          coord = []
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

  return {'species':species, 'lattice':np.array(cell_params), 'abc':np.array(abc)}


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
      block = mf.pop(0)
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

  if struct['aunit'] == 'angstrom':
    struct['abc'] = struct['abc'] @ np.linalg.inv(struct['lattice'])

  return struct


def struct_from_inputfile_CIF ( fname:str ):
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
        #if '_cell_len' in tag:
        #  ind = ord(tag.split('_')[-1]) - 97
        #  #cell_param[ind] = strip_float(lines[eL+1])
        #  data[tag] = strip_float(lines[eL+1])
      tag = lambda t : t.split('_')[-1]
      if '_cell_length_' in lines[eL]:
        ls = lines[eL].split()
        ind = ord(tag(ls[0])) - 97
        cell_param[ind] = strip_float(ls[1])
      elif '_cell_angle_' in lines[eL]:
        ls = lines[eL].split()
        ind = 3 + ['alpha', 'beta', 'gamma'].index(tag(ls[0]))
        cell_param[ind] = strip_float(ls[1])
      else:
        print(tag(lines[eL].split()[0]))
      eL += 1
      #data_blocks.append({'id':prefix, 'cell_param':cell_param, '1':1})
      #cell_param = np.zeros(6, dtype=float)
      '''

  print(data_blocks)
  quit()
    



def infer_file_type ( fname:str ):

  extension = fname.split('.')
  if len(extension) == 1:
    print('Cannot infer file type without an extension. Assuming qe')
    return 'qe'
  extension = extension[-1].lower()
  if extension == 'in' or extension == 'out':
    ftype = 'qe'
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
    else:
      raise ValueError('Cannot read file type {}'.format(ftype))

  except Exception as e:
    print('Failed to read file {}'.format(fname))
    raise e

