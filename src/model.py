
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
    

