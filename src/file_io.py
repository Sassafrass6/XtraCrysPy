from numpy import ndarray


def struct_from_outputfile_QE ( fname:str ):
  '''
  '''
  from os.path import isfile,join
  import numpy as np

  if not isfile(fname):
    msg = 'File {} does not exist.'.format(join(getcwd(),fname))
    raise FileNotFoundError(msg)

  struct = {'lunit':'bohr', 'aunit':'alat'}
  with open(fname, 'r') as f:
    lines = f.readlines()

    eL = 0
    nL = len(lines)

    try:
      struct['species'] = []
      celldm = np.empty(6, dtype=float)
      while 'bravais-lattice' not in lines[eL]:
        eL += 1
      ibrav = int(lines[eL].split()[3])

      while  'celldm' not in lines[eL]:
        eL += 1
      celldm[:3] = [float(v) for i,v in enumerate(lines[eL].split()) if i%2==1]
      celldm[3:] = [float(v) for i,v in enumerate(lines[eL+1].split()) if i%2==1]

      if ibrav != 0:
        from .lattice_format import lattice_format_QE
        struct['lattice'] = lattice_format_QE(ibrav, celldm)
      else:
        while 'crystal axes' not in lines[eL]:
          eL += 1
        coord = []
        for l in lines[eL+1:eL+4]:
          coord.append([celldm[0]*float(v) for v in l.split()[3:6]])
        struct['lattice'] = np.array(coord)

      while 'site n.' not in lines[eL]:
        eL += 1
      eL += 1
      apos = []
      while 'End' not in lines[eL] and lines[eL] != '\n':
        line = lines[eL].split()
        struct['species'].append(line[1])
        apos.append([float(v) for v in line[6:9]])
        eL += 1
      apos = celldm[0] * np.array(apos)
      struct['abc'] = apos @ np.linalg.inv(struct['lattice'])

    except Exception as e:
      print('ERROR: Could not read the QE output.')
      raise e

  return struct


def read_relaxed_coordinates_QE ( fname:str, read_all:bool=True ):
  '''
    Reads relaxed atomic positions from a QE .out file. If CELL_PARAMETERS is present, the crystal coordinates are also read.

    Arguments:
      fname (str): File name (including path) for the .out file
      read_all (bool): True forces all relax steps to be read. If EoF is encountered before 'final coordinates' the last coordinates to appear in the file are retunred. If no coordinates are found, an empty dictionary is returned.

    Returns:
      (dict): Dictionary with one or two entries - 'apos' for atomic positions and 'coord' for crystal coordinates.
  '''
  from .conversion import ANG_BOHR
  from os.path import isfile,join
  from os import getcwd
  import numpy as np
  import re

  abc = []
  cell_params = []
  struct = struct_from_outputfile_QE(fname)

  with open(fname, 'r') as f:
    lines = f.readlines()

    eL = 0
    nL = len(lines)

    try:
      def read_apos ( sind ):
        unit = lines[sind].split()[1].strip('(){{}}')
        apos = []
        sind += 1
        while lines[sind] != '\n' and not 'End final coordinates' in lines[sind]:
          apos.append([float(v) for v in lines[sind].split()[1:4]])
          sind += 1
        apos = np.array(apos)

        if len(unit) > 1:
          struct['aunit'] = unit
        if 'angstrom' in unit:
          apos = ANG_BOHR * (apos @ np.linalg.inv(struct['lattice']))
        elif 'bohr' in unit:
          apos = apos @ np.linalg.inv(struct['lattice'])
        return sind, apos

      while eL < nL:
        while eL < nL and 'CELL_PARAMETERS' not in lines[eL] and 'ATOMIC_POSITIONS' not in lines[eL]:
          eL += 1
        if eL >= nL:
          break

        if 'ATOMIC_POSITIONS' in lines[eL]:
          sind,apos = read_apos(eL)
          abc.append(apos)
          eL = sind

        elif ('CELL_PARAMETERS' in lines[eL]):
          coord = []
          unit = lines[eL].split()[1].strip('(){{}}')
          alat = 1
          if 'alat' in unit or len(unit) == 0:
            struct['lunit'] = 'alat'
            if 'alat' in unit:
              cpattern = re.search('\(([^\)]+)\)', lines[eL])
              if cpattern is not None:
                alat = float(cpattern.group(0)[1:-1].split('=')[1])
          else:
            struct['lunit'] = unit
          for l in lines[eL+1:eL+4]:
            coord.append(alat*np.array([float(v) for v in l.split()]))
          cell_params.append(coord)
          eL += 4

          while 'ATOMIC_POSITIONS' not in lines[eL]:
            eL += 1
          eL,apos = read_apos(eL)
          abc.append(apos)

    except Exception as e:
      print('WARNING: No atomic positions or cell coordinates were found.', flush=True)
      raise e

  # If it's a vc-relax keep only the last vectors+positions
  # If it's a relax, keep only the last positions, cell vectors from header
  # else, pass the initial + full relaxation trajectory
  if len(cell_params)>0 and ( not read_all ):
    cell_params = np.array(cell_params[-1])
    abc = abc[-1]
  elif ( not read_all ):
    cell_params = struct['lattice']
    abc = abc[-1]
  else:
    cell_params = np.array([struct['lattice']] + cell_params)
    abc = np.array([struct['abc']] + abc)

  struct['lunit'] = 'bohr'
  struct['aunit'] = 'crystal'
  struct['lattice'] = cell_params
  struct['abc'] = abc

  return struct


def read_relaxed_coordinates_CP2K_XYZ ( fname:str ):
  import numpy as np

  lines = None
  with open(fname, 'r') as f:
    lines = f.readlines()

  i = 0
  abc = []
  enes = []
  specs = []
  nl = len(lines)
  while i < nl:
    nat = int(lines[i].strip())
    ls = lines[i+1].split()
    niter = int(ls[2].strip(','))
    enes.append(float(ls[-1]))
    i += 2
    pos = []
    for j in range(nat):
      ls = lines[i+j].split()
      if i < nat:
        specs.append(ls[0])
      pos.append([float(v) for v in ls[1:]])
    i += nat
    abc.append(pos)

  struct = {}
  struct['species'] = specs
  struct['abc'] = 1.88973 * np.array(abc)
  struct['lattice'] = np.array([np.eye(3)], dtype=float)
  for i in range(3):
    struct['lattice'][0,i] *= np.max(struct['abc'][:,:,i])
  linv = np.linalg.inv(struct['lattice'])
  for i,p in enumerate(struct['abc']):
    struct['abc'][i] = p @ linv

  return struct


def md_coordinates_LAMMPS ( fname:str ):
  '''
  '''
  import numpy as np

  lines = None
  with open(fname, 'r') as f:
    lines = f.readlines()

  il = 0
  ll = len(lines)
  species = None
  positions = []
  while il < ll:
    nat = int(lines[il].split()[0])
    il += 2
    atoms = []
    species = []
    for _ in range(nat):
      ls = lines[il].split()
      spec = int(ls[0]) if ls[0].isnumeric() else ls[0]
      species.append(spec)
      atoms.append([float(v) for v in ls[1:]])
      il += 1
    positions.append(atoms)

  struct = {}
  struct['species'] = species
  struct['abc'] = np.array(positions)
  struct['lattice'] = np.array([np.eye(3)], dtype=float)
  for i in range(3):
    struct['lattice'][0,i] *= np.max(struct['abc'][:,:,i])
  linv = np.linalg.inv(struct['lattice'])
  for i,p in enumerate(struct['abc']):
    struct['abc'][i] = p @ linv

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
  from .conversion import ANG_BOHR
  from os.path import isfile
  import numpy as np
  import re

  if not isfile(fname):
    raise FileNotFoundError('File {} does not exist.'.format(fname))

  fstr = None
  with open(fname, 'r') as f:
    fstr = f.read()

  # Datatype format helpers for QE input
  nocomma = lambda s : s.replace(',', '')
  qebool = lambda s : True if s.split('.')[1][0].lower() == 't' else False
  qenum = lambda s : s.split('=')[1].replace('d', 'e')
  qeint = lambda s : int(qenum(s))
  qefloat = lambda s : float(qenum(s))
  def inquote ( s ):
    v = '"' if '"' in s else "'"
    return s.split(v)[1]

  # Process blocks
  struct = {}
  pattern = re.compile('&(.*?)/@')
  celldm = np.zeros(6, dtype=float)
  comment = lambda v : v != '' and '!' not in v
  matches = pattern.findall(fstr.replace(' ', '').replace('\n', '@ '))
  for m in matches:
    m = [s.replace(' ', '').split('!')[0] for s in re.split(', |@', m)]
    mf = []
    for v in m:
      mf += v.split(',')
    block = mf.pop(0).lower()
    mf = filter(comment, mf)
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
        elif 'nat' in s:
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
        else:
          cpat = ['A=', 'B=', 'C=', 'cosAB', 'cosAC', 'cosBC']
          for i,p in enumerate(cpat):
            if p in s:
              celldm[i] = ANG_BOHR * qefloat(s)
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

  # Process CARDS
  fstr = list(filter(comment, fstr.split('\n')))
  def scan_blank_lines ( nl ):
    nl += 1
    while fstr[nl] == '':
      nl += 1
    return nl

  il = 0
  cl = 0
  kl = 0
  while 'ATOMIC_POSITION' not in fstr[il]:
    if 'CELL_PARAM' in fstr[il]:
      cl = il
    if 'K_POINTS' in fstr[il]:
      kl = il
    il += 1
  unit = fstr[il].split()
  struct['aunit'] = unit[1].strip('(){{}}').lower() if len(unit) > 1 else 'alat'

  il = scan_blank_lines(il)

  spec = []
  abc = np.empty((struct['nat'],3), dtype=float)
  for i in range(struct['nat']):
    ls = fstr[il+i].split()
    spec.append(ls[0])
    abc[i,:] = np.array([float(v) for v in ls[1:4]])

  struct['abc'] = abc
  struct['species'] = spec

  if kl == 0:
    while 'K_POINTS' not in fstr[kl]:
      if 'CELL_PARAM' in fstr[cl]:
        cl = kl
      kl += 1

  if 'gamma' in fstr[kl].lower():
    struct['kpnts'] = {'type':'gamma'}
  elif 'automatic' in fstr[kl]:
    kl = scan_blank_lines(kl)
    nk1,nk2,nk3,o1 = (int(v) for v in fstr[kl].split()[:4])
    struct['kpnts'] = {'type':'automatic', 'nk1':nk1, 'nk2':nk2, 'nk3':nk3, 'offset':o1}

  nf = len(fstr)
  if cl == 0:
    while cl < nf and 'CELL_PARAM' not in fstr[cl]:
      cl += 1

  if cl < nf:
    unit = fstr[cl].split()
    struct['lunit'] = unit[1].strip('(){{}}').lower() if len(unit) > 1 else 'alat'
    cl = scan_blank_lines(cl)
    struct['lattice'] = np.array([np.array([float(v) for v in fstr[cl+c].split()]) for c in range(3)])

  if 'lattice' not in struct:
    from .lattice_format import lattice_format_QE
    struct['lattice'] = lattice_format_QE(struct['ibrav'], struct['celldm'])

  if struct['lunit'] == 'angstrom':
    struct['lattice'] *= ANG_BOHR
  elif struct['lunit'] == 'alat':
    struct['lattice'] *= celldm[0]
  elif struct['lunit'] != 'bohr':
    print('WARNING: Unit {} may behave strangely'.format(struct['lunit']))
  struct['lunit'] = 'bohr'

  if struct['aunit'] in ['bohr', 'angstrom', 'alat']:
    if struct['aunit'] == 'angstrom':
      struct['abc'] *= ANG_BOHR
    elif struct['aunit'] == 'alat':
      struct['abc'] = struct['abc'] @ struct['lattice']
    struct['abc'] = struct['abc'] @ np.linalg.inv(struct['lattice'])
  struct['aunit'] = 'crystal'

  return struct


def struct_from_atoms_ASE ( atoms ):
  from .conversion import ANG_BOHR
  import numpy as np

  syms = atoms.symbols
  la = len(syms)

  def extract_number ( ind ):
    nn = 1
    while ind+nn < la and syms[ind+nn].isnumeric():
      nn += 1
    return nn-1, int(syms[ind:ind+nn])

  ci = 0
  spec = []
  while ci < la:
    if ci == la-1:
      spec.append(syms[ci])
      ci += 1
    else:
      num = 1
      sym = syms[ci]
      if syms[ci+1].isnumeric():
        ci += 1
        nn,num = extract_number(ci)
        ci += nn
      elif syms[ci+1].islower():
        ci += 2
        sym += syms[ci-1]
        if ci < la and syms[ci].isnumeric():
          nn,num = extract_number(ci)
          ci += nn
      else:
        ci += 1
      spec += num * [sym]

  if np.all(atoms.cell == 0.0):
    msg = 'ASE failed to read atomic information.'
    raise Exception(msg)

  lat = np.empty((3,3), dtype=float)
  for i,v in enumerate(atoms.cell):
    lat[i,:] = v[:]

  abc = atoms.positions.copy() @ np.linalg.inv(lat)
  lat *= ANG_BOHR

  return lat, spec, abc


def struct_from_inputfile_ASE ( fname:str, format=None, index=None ):
  '''
  '''
  import ase.io as aio
  import numpy as np

  atoms = aio.read(fname, format=format, index=index)

  if isinstance(atoms, list):
    if len(atoms) == 0:
      raise Exception(f'No structures were read from file {fname}.')

    species = None
    lattice,positions = [],[]
    for a in atoms:
      lat,spec,abc = struct_from_atoms_ASE(a)
      if species is None:
        species = spec
      lattice.append(lat)
      positions.append(abc)

    return {'lattice':np.array(lattice), 'species':species, 'abc':np.array(positions)}

  else:
    lat,spec,abc = struct_from_atoms_ASE(atoms)

    return {'lattice':lat, 'species':spec, 'abc':abc}


def struct_from_inputfile_CP2K ( fname:str ):
  '''
  '''
  from .conversion import ANG_BOHR
  import numpy as np
  import re

  lines = None
  with open(fname, 'r') as f:
    lines = f.readlines()

  isl = 0
  nl = len(lines)
  while '&SUBSYS' not in lines[isl]:
    isl += 1
  isl += 1
  ssl = isl

  def unit_conv ( unit, angular=False ):
    if angular:
      if unit is None or 'deg' in unit.group(0).lower():
        return np.pi / 180
    else:
      if unit is None:
        return ANG_BOHR
      else:
        unit = unit.group(0).lower()
        if 'bohr' not in unit:
          if 'angstrom' in unit:
            return ANG_BOHR
          elif 'pm' in unit:
            return 100 * ANG_BOHR
          elif 'nm' in unit:
            return 0.1 * ANG_BOHR
          else:
            return 1e-10 * ANG_BOHR
    return 1

  struct = {}
  crystal_positions = False
  while 'SUBSYS' not in lines[isl]:

    if '&CELL' in lines[isl]:
      nc = 1
      cell = {}
      lattice = np.zeros((3,3), dtype=float)
      while 'CELL' not in lines[nc+isl]:
        unit = re.search('\[([^.]+)\]', lines[nc+isl])
        angular = 'ALPHA' in lines[nc+isl]
        conv = unit_conv(unit, angular=angular)
        if unit is not None:
          lines[nc+isl] = lines[nc+isl].replace(unit.group(0), '')
        ls = lines[nc+isl].split()
        if len(ls) > 2:
          cell[ls[0]] = conv * np.array([float(v) for v in ls[1:]])
        else:
          cell[ls[0]] = ls[1].lower() if len(ls)==2 else 'cp2k'
        nc += 1
      isl += nc

      abg = None
      abc = ['A', 'B', 'C']
      cff = 'CELL_FILE_FORMAT'
      if cff in cell and not cell[cff] == 'cp2k':
        lfname = cell['CELL_FILE_NAME']
        lattice = struct_from_inputfile_ASE(lfname)['lattice']
      elif 'ALPHA_BETA_GAMMA' in cell:
        from .lattice_format import lattice_format_abc_abg
        if 'ABC' in cell:
          cABC = cell['ABC']
        else:
          cABC = [cell[c] for c in abc]
        abg = cell['ALPHA_BETA_GAMMA']
        args = [cell[c] for c in cABC] + [v for v in abg]
        lattice[:,:] = lattice_format_abc_abg(*args)
      else:
        for i,c in enumerate(abc):
          lattice[i,:] = cell[c]
      struct['lattice'] = lattice

    elif '&COORD' in lines[isl]:
      nc = 2
      conv = 1
      species = []
      positions = []
      ls = lines[1+isl].split()
      if 'SCALED' in ls:
        crystal_positions = True
        if len(ls) > 1 and 'false' in ls[1].lower():
          raise Exception('Specify units with UNIT tag')
      elif 'UNIT' in ls:
        if len(ls) > 1:
          conv = unit_conv(ls[1])
      while 'COORD' not in lines[nc+isl]:
        ls = lines[nc+isl].split()
        species.append(ls[0])
        positions.append([float(v) for v in ls[1:4]])
        nc += 1
      positions = conv * np.array(positions)
      struct['abc'] = positions
      struct['species'] = species
      isl += nc

    else:
      isl += 1

  if not crystal_positions:
    struct['abc'] = struct['abc'] @ np.linalg.inv(struct['lattice'])

  return struct


def struct_from_inputfile ( fname:str, ftype=None, index=None ):
  '''
  '''

  try:
    if ftype == 'cp2k-in':
      return struct_from_inputfile_CP2K(fname)

    elif ftype == 'cp2k-xyz':
      return read_relaxed_coordinates_CP2K_XYZ(fname)

    elif ftype == 'lammps-traj':
      return md_coordinates_LAMMPS(fname)

    else:
      try:
        return struct_from_inputfile_ASE(fname, format=ftype, index=index)
      except Exception as e:
        print('Attempting read with internal QE parser.\n')
        if ftype == 'espresso-out':
          return read_relaxed_coordinates_QE(fname)
        else:
          return struct_from_inputfile_QE(fname)

  except Exception as e:
    print('Failed to read file {}'.format(fname))
    raise e


def struct_from_file_sequence ( fnames ):
  '''
  '''
  import numpy as np

  nfn = len(fnames)
  if nfn == 0:
    raise Exception('File list is empty')

  struct = struct_from_inputfile(fnames[0])

  spec = struct['species']
  lats = np.empty((nfn,3,3), dtype=float)
  abcs = np.empty((nfn,len(spec),3), dtype=float)
  lats[0] = struct['lattice']
  abcs[0] = struct['abc']
  for i in range(1,nfn):
    struct = struct_from_inputfile(fnames[i])
    lats[i] = struct['lattice']
    abcs[i] = struct['abc']

  struct['lattice'] = lats
  struct['abc'] = abcs

  return struct


def read_CUBE ( fname:str ):
  '''
  '''
  import numpy as np

  lines = None
  with open (fname, 'r') as f:
    lines = f.readlines()

  pop = lambda : lines.pop(0)
  while not lines[0].split()[0].isnumeric():
    pop()

  pop()
  grid = np.empty(3, dtype=int)
  for i in range(3):
    grid[i] = int(pop().split()[0])
  for _ in range(3):
    pop()

  data = np.empty(grid, dtype=float)
  for i in range(grid[0]):
    for j in range(grid[1]):
      lcnt = 0
      while lcnt < grid[2]:
        tdata = [float(v) for v in pop().split()]
        data[i,j,lcnt:lcnt+len(tdata)] = np.array(tdata)
        lcnt += len(tdata)
  return data[::-1, ::-1, ::-1]


def read_XSF ( fname:str ):
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
        dims = np.array([int(v) for v in lines[ind+1].split()])
        ind += 6
        data = np.empty(dims, dtype=float)
        for i in range(dims[0]):
          for j in range(dims[1]):
            data[i,j,:] = np.array([float(v) for v in lines[ind].split()])
            ind += 1
        blocks.append(data)
      ind += 1

  return np.array(blocks)
