from fury.utils import colors_from_actor, update_actor
from fury.utils import vertices_from_actor
from .Atomic import Atomic
import numpy as np

from time import time

class MolDyn ( Atomic ):

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True,
                 background=(0,0,0), perspective=False, model=None,
                 params={}, multi_frame=False, nsc=(1,1,1),
                 bond_type='Stick', sel_type='Chain', unit='angstrom',
		 runit='degree', constrain_atoms=False, atom_res=(0,0),
                 image_prefix='XCP_Image', resolution=4, dt=0.1):
    super().__init__(size, axes, boundary, background,
                     perspective, model, params, multi_frame, nsc, bond_type, sel_type, unit, runit, constrain_atoms, atom_res, image_prefix, resolution, 0)

    self.fn = 0
    self.dt = dt
    self.velocities = np.zeros((self.model.atoms.shape[0],3), dtype=float)
    self.st = time()
    self.smanager.add_timer_callback(True, 30, self.frame_timer)


  def lj_forces ( self, sigma=0.8, epsilon=.05 ):
    atoms = self.model.atoms
    species = self.model.species

    forces = np.zeros((atoms.shape[0],3), dtype=float)
    for i in range(atoms.shape[0]-1):
      for j in range(i+1, atoms.shape[0]):
        s1,s2 = species[i],species[j]
        r = atoms[j] - atoms[i]
        rm = np.linalg.norm(r)
        fij = 24*epsilon/sigma * (2*(sigma/rm)**13-(sigma/rm)**7)
        fij *= r/rm
        forces[j] += fij
        forces[i] -= fij

    return forces 


  def update_atoms ( self ):
    forces = self.lj_forces()
    atoms = self.model.atoms
    for i in range(atoms.shape[0]):
      self.velocities[i] += forces[i] * self.dt
      atoms[i] += self.velocities[i] * self.dt
    self.update_atomic_positions()

  def frame_timer ( self, obj, event ):
    print(f'Frame {self.fn} at time ', time() - self.st)
    self.update_atoms()
    self.fn += 1


