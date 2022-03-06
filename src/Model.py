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
      self.units = 'bohr'
      self.bond_type = None
      self.atoms = params['abc']
      self.species = params['species']
      self.lattice = params['lattice']
    except Exception as e:
      print('Manual structures require \'lattice\', \'species\', and atomic positions \'abc\'')
      raise e

    if 'units' in params:
      self.units = params['units']

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

    if 'units' in params:
      unit = params['units']
      if unit == 'bohr':
        pass
      elif unit == 'angstrom':
        from .conversion import ANG_BOHR as conv
        self.lattice *= conv
        for k in self.bonds.keys():
          self.bonds[k] *= conv

    mins = np.min(self.atoms, axis=0)
    maxs = np.max(self.atoms, axis=0)
    if not np.all([mi>=-1 and ma<=1 for mi,ma in zip(mins,maxs)]):
      print('WARNING: \'abc\' requires crystal coordinates. Specify \'xyz\' for custom units..')

    for s in self.species:
      if not s in self.colors:
        if s in atom_defaults:
          self.colors[s] = np.array(atom_defaults[s]['color'])/255
        else:
          self.colors[s] = (.8,.8,.8)
      if s not in self.radii:
        if s in atom_defaults:
          self.radii[s] = atom_defaults[s]['radius']
        else:
          self.radii[s] = 1


  def lattice_atoms_bonds ( self, nc1, nc2, nc3, bond_type='stick' ):
    '''
    '''

    nsc = (nc1,nc2,nc3)
    nat = self.atoms.shape[0]

    if self.bond_type is None:
      self.bond_type = bond_type

    lattice = self.lattice.copy()
    for i,n in enumerate(nsc):
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
    for i,n in enumerate(nsc):
      atoms[:,i] /= n
    atoms = atoms @ lattice

    bonds = []
    bdirs = []
    bcols = []
    brads = []
    bheight = []
    natsc = atoms.shape[0]
    skey = lambda v1,v2 : v1 + '_' + v2
    if bond_type != 'sphere':
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
                brad = .5 / dist
              elif bond_type == 'primary':
                brad = 1.5 * prad / dist
              else:
                raise ValueError('Bond types are \'stick\' and \'primary\'')
              brads.append(brad)
              bheight.append(dist/2)
    bonds = np.array(bonds)
    bdirs = np.array(bdirs)
    bcols = np.array(bcols)
    brads = np.array(brads)
    bheight = np.array(bheight)

    lines = []
    lat = lattice
    for ix in range(nc1):
      for iy in range(nc2):
        for iz in range(nc3):
          orig = np.array([ix,iy,iz]) @ lat
          corner = np.sum(lat, axis=0) + orig
          lines += [[orig,orig+a] for a in lat]
          lines += [[orig+p,corner] for p in [lat[i]+lat[j] for i in range(3) for j in range(i+1,3)]]
          lines += [[orig+lat[k],orig+lat[i]+lat[j]] for i in range(3) for j in range(i+1,3) for k in [i,j]]
    lines = np.array(lines)
    for i,n in enumerate(nsc):
      lines[:,:,i] /= n

    return (atoms,acols,radii), (bonds,bdirs,bcols,brads,bheight), (lines)
