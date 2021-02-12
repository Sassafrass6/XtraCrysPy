
def get_bond_pairs ( natoms, atoms, bonds, labels ):
  '''
  '''
  import numpy as np

  bond_pairs = []
  if bonds is None or len(atoms) < 2:
    return bond_pairs

  na = natoms
  for i,k in enumerate(bonds.keys()):
    sA,sB = tuple(k.split('_'))
    for m in range(len(atoms)):
      for n in range(m+1, len(atoms)):
        l1,l2 = labels[m%na],labels[n%na]
        if (l1 == sA and l2 == sB) or (l1 == sB and l2 == sA):
          if np.linalg.norm(atoms[m] - atoms[n]) <= bonds[k]:
            bond_pairs += [[m, n]]
  return bond_pairs

cube_species_mapping = {
  1: 'H',
  2: 'He',
  3: 'Li',
  4: 'Be',
  5: 'B',
  6: 'C',
  7: 'N',
  8: 'O',
  9: 'F',
  10: 'Ne',
  11: 'Na',
  12: 'Mg',
  13: 'Al',
  14: 'Si',
  15: 'P',
  16: 'S',
  17: 'Cl',
  18: 'Ar',
  19: 'K',
  20: 'Ca',
  21: 'Sc',
  22: 'Ti',
  23: 'V',
  24: 'Cr',
  25: 'Mn',
  26: 'Fe',
  27: 'Co',
  28: 'Ni',
  29: 'Cu',
  30: 'Zn',
  31: 'Ga',
  32: 'Ge',
  33: 'As',
  34: 'Se',
  35: 'Br',
  36: 'Kr',
  37: 'Rb',
  38: 'Sr',
  39: 'Y',
  40: 'Zr',
  41: 'Nb',
  42: 'Mo',
  43: 'Tc',
  44: 'Ru',
  45: 'Rh',
  46: 'Pd',
  47: 'Ag',
  48: 'Cd',
  49: 'In',
  50: 'Sn',
  51: 'Sb',
  52: 'Te',
  53: 'I',
  54: 'Xe',
  55: 'Cs',
  56: 'Ba',
  57: 'La',
  58: 'Ce',
  59: 'Pr',
  60: 'Nd',
  61: 'Pm',
  62: 'Sm',
  63: 'Eu',
  64: 'Gd',
  65: 'Tb',
  66: 'Dy',
  67: 'Ho',
  68: 'Er',
  69: 'Tm',
  70: 'Yb',
  71: 'Lu',
  72: 'Hf',
  73: 'Ta',
  74: 'W',
  75: 'Re',
  76: 'Os',
  77: 'Ir',
  78: 'Pt',
  79: 'Au',
  80: 'Hg',
  81: 'Tl',
  82: 'Pb',
  83: 'Bi',
  84: 'Po',
  85: 'At',
  86: 'Rn',
  87: 'Fr',
  88: 'Ra',
  89: 'Ac',
  90: 'Th',
  91: 'Pa',
  92: 'U',
  93: 'Np',
  94: 'Pu',
  95: 'Am',
  96: 'Cm',
  97: 'Bk',
  98: 'Cf',
  99: 'Es',
  100: 'Fm',
  101: 'Md',
  102: 'No',
  103: 'Lr',
  104: 'Rf',
  105: 'Db',
  106: 'Sg',
  107: 'Bh',
  108: 'Hs',
  109: 'Mt',
  110: 'Ds',
  111: 'Rg',
  112: 'Cn',
  113: 'Nh',
  114: 'Fl',
  115: 'Mc',
  116: 'Lv',
  117: 'Ts',
  118: 'Og',
}

# radii in angstrom
cpk_species = {
  'C': {
    'color': (0.784, 0.784, 0.784, 0),
    'radius': 0.963
  },
  'O': {
    'color': (0.941, 0, 0, 0),
    'radius': 0.875
  },
  'H': {
    'color': (1, 1, 1, 0),
    'radius': 0.722
  },
  'N': {
    'color': (0.561, 0.561, 1, 0),
    'radius': 0.917
  },
  'S': {
    'color': (1, 0.784, 0.216, 0),
    'radius': 1.009
  },
  'P': {
    'color': (1, 0.647, 0, 0),
    'radius': 1.037
  },
  'Cl': {
    'color': (0, 1, 0, 0),
    'radius': 0.989
  },
  'Br': {
    'color': (0.647, 0.165, 0.165, 0),
    'radius': 1.047
  },
  'Zn': {
    'color': (0.647, 0.165, 0.165, 0),
    'radius': 0.691
  },
  'Ni': {
    'color': (0.647, 0.165, 0.165, 0),
    'radius': 0.708
  },
  'Cu': {
    'color': (0.647, 0.165, 0.165, 0),
    'radius': 0.874
  },
  'Na': {
    'color': (0, 0, 1, 0),
    'radius': 0.841
  },
  'Fe': {
    'color': (1, 0.647, 0, 0),
    'radius': 0.728
  },
  'Mg': {
    'color': (0.164, 0.5, 0.164, 0),
    'radius': 0.755
  },
  'Ca': {
    'color': (0.5, 0.5, 0.5, 0),
    'radius': 0.850
  },
}

organic_bonds = {
  'C_H': 1.12*2,
  'C_C': 1.54*2,
  'C_N': 2.10*2,
  'C_O': 2.15*2,
  'N_H': 1.10*2,
  'N_N': 1.55*2,
  'O_H': 1.00*2,
  'O_O': 1.54*2
}