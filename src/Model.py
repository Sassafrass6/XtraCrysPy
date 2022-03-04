import numpy as np

class Model:

  def __init__ ( self, params={}, fname=None, relax=False ):
    '''

    Arguments:
    '''
    from .defaults import atom_defaults

    if fname is not None:
      if relax:
        # Read DFT output file
        pass
      else:
        from .dft_file_io import struct_from_inputfile
        params.update(struct_from_inputfile(fname))

    try:
      self.atoms = params['abc']
      self.species = params['species']
      self.lattice = params['lattice']
    except Exception as e:
      print('Manual structures require \'lattice\', \'species\', and atomic positions \'abc\'')
      raise e

    self.bonds = {}
    if 'bonds' in params:
      self.bonds = params['bonds']

    self.colors = {}
    if 'colors' in params:
      for k,v in params['colors'].items():
        self.colors[k] = np.array(v)

    self.radii = {}
    if 'radii' in params:
      for k,v in params['radii'].items():
        self.radii[k] = v

    for s in self.species:
      if not s in self.colors:
        if s in atom_defaults:
          self.radii[s] = atom_defaults[s]['radius']
          self.colors[s] = atom_defaults[s]['color']
        else:
          self.radii[s] = 1
          self.colors[s] = (.8,.8,.8)


  def lattice_atoms_bonds ( self, nc1, nc2, nc3, bond_type='stick' ):
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
            spec = self.species[i]
            atoms.append(origin + a)
            radii.append(self.radii[spec])
            acols.append(self.colors[spec])
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
    if bond_type == 'primary':
      prad = np.min([r for k,r in self.radii.items()])
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
            if bond_type == 'stick':
              brad = .75 / dist
            elif bond_type == 'primary':
              brad = 2 * prad / dist
            else:
              raise ValueError('Bond types are \'stick\' and \'primary\'')
            brads.append(brad)
            bheight.append(dist/2)
    bonds = np.array(bonds)
    bdirs = np.array(bdirs)
    bcols = np.array(bcols)
    brads = np.array(brads)
    bheight = np.array(bheight)

    return (atoms,acols,radii), (bonds,bdirs,bcols,brads,bheight)
