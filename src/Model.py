import numpy as np

class Model:

  def __init__ ( self, structure ):
    '''

    Arguments:
    '''

    self.atoms = structure['abc']
    self.species = structure['species']
    self.lattice = structure['lattice']

    self.bonds = {}
    if 'bonds' in structure:
      self.bonds = structure['bonds']

    self.colors = {}
    if 'colors' in structure:
      for k,v in structure['colors'].items():
        self.colors[k] = np.array(v)
    # AUTO FILL DEFAULTS

    self.radii = {}
    if 'radii' in structure:
      for k,v in structure['radii'].items():
        self.radii[k] = v
    # AUTO FILL DEFAULTS


  def lattice_atoms_bonds ( self, nc1, nc2, nc3 ):
    '''
    '''

    nat = self.atoms.shape[0]

    lattice = self.lattice.copy()
    for i,n in enumerate([nc1,nc2,nc3]):
      lattice[:,i] *= n

    atoms = []
    acols = []
    radii = []
    for n1 in range(nc1):
      for n2 in range(nc2):
        for n3 in range(nc3):
          origin = np.array([n1,n2,n3])
          for i,a in enumerate(self.atoms):
            atoms.append(origin + a)
            acols.append(self.colors[self.species[i]])
            radii.append(1)
    atoms = np.array(atoms)
    acols = np.array(acols)
    radii = np.array(radii)
    for i,n in enumerate([nc1,nc2,nc3]):
      atoms[:,i] /= n
    atoms = atoms @ lattice

    bonds = []
    bdirs = []
    bcols = []
    brads = []
    bheight = []
    natsc = atoms.shape[0]
    skey = lambda v1,v2 : v1 + '_' + v2
    for i1 in range(natsc-1):
      for i2 in range(i1+1, natsc):
        a1 = atoms[i1]
        a2 = atoms[i2]
        t1 = self.species[i1%nat]
        t2 = self.species[i2%nat]

        bnd_len = None
        if skey(t1,t2) in self.bonds:
          bnd_len = self.bonds[skey(t1,t2)]
        elif skey(t2,t1) in self.bonds:
          bnd_len = self.bonds[skey(t2,t1)]
        
        if not bnd_len is None:
          conn = a2 - a1
          dist = np.linalg.norm(conn)
          if bnd_len >= dist:
            cent = (a1 + a2) / 2
            bonds.append([cent-conn/4, cent+conn/4])
            bdir = conn/dist
            bdirs.append([bdir]*2)
            bcols.append([self.colors[t1], self.colors[t2]])
            brads.append((1+1)/8)
            bheight.append(dist/2)
    bonds = np.array(bonds)
    bdirs = np.array(bdirs)
    bcols = np.array(bcols)
    brads = np.array(brads)
    bheight = np.array(bheight)

    return (atoms,acols,radii), (bonds,bdirs,bcols,brads,bheight)
