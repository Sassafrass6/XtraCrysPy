from numpy import ndarray


def read_relaxed_coordinates ( fname:str, vcrelax=False, read_all=False ):
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

  apos = []
  coord = [] if vcrelax else None 
  with open(fname, 'r') as f:
    lines = f.readlines()

    eL = 0
    nL = len(lines)

    if not read_all:
      while eL < nL and 'Begin final coordinates' not in lines[eL]:
        eL += 1
      if eL == nL:
        raise EOFError('Reached end of file before finding final coordinates.')

    try:
      while eL < nL:
        if vcrelax:
          while 'CELL_PARAMETERS' not in lines[eL]:
            eL += 1
          coord = []
          for l in lines[eL+1:eL+4]:
            coord.append(np.array([float(v) for v in l.split()]))
          eL += 4

        while 'ATOMIC_POSITIONS' not in lines[eL]:
          eL += 1
        eL += 1

        apos = []
        while lines[eL] != '\n' and not 'End final coordinates' in lines[eL]:
          apos.append(np.array([float(v) for v in lines[eL].split()[1:4]]))
          eL += 1

        if 'End final coordinates' in lines[eL]:
          break
    except Exception as e:
      if not vcrelax:
        print('WARNING: Reached EoF without finding Final Coordinates.', flush=True)
      if len(apos) == 0:
        print('WARNING: No atomic positions or cell coordinates were found.', flush=True)
        return {}

  cdict = {'apos':np.array(apos)}
  if vcrelax:
    cdict['coord'] = np.array(coord)

  return cdict


def struct_from_inputfile ( fname:str ) -> dict:
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

