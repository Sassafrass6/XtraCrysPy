from .XtraCrysPy import XtraCrysPy
from .Model import Model
import numpy as np

class XCP_BZ ( XtraCrysPy ):

  def __init__ ( self, size=(1024, 1024), axes=True, perspective=False, model=None, frame_width=4e-3 ):
    super().__init__(size, axes, perspective)

    self.model = model
    self.frame_width = frame_width
    if model is not None:
      if isinstance(model, dict):
        if 'rlattice' in model:
          self.rlattice = model['rlattice']
        elif 'lattice' in model:
          print('Converting lattice to rlattice')
          lat = model['lattice']
          self.rlattice = np.empty((3,3), dtype=float)
          for i in range(3):
            self.rlattice[i,:] = np.cross(lat[i-2], lat[i-1])
          volume = lat[0].dot(np.cross(lat[1], lat[2]))
          self.rlattice *= 2 * np.pi / volume
        else:
          raise KeyError('No reciprocal lattice (rlattice) or lattice (lattice) specified in model dictionary.')

      elif isinstance(model, Model):
        self.rlattice = model.rlattice

      else:
        raise ValueError('Argument \'model\' must be a Model object or a dictionary.')

      self.bravais_boundaries()


  def bravais_boundaries ( self ):
    from numpy.linalg import det,norm,solve
    from fury import actor

    b_vec = self.rlattice
    indices = [[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],[-1,0,-1],[1,1,0],[-1,-1,0],[1,1,1],[-1,-1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,-1],[-1,1,1],[1,-1,-1],[-1,1,0],[1,-1,0],[1,0,-1],[-1,0,1],[0,1,-1],[0,-1,1]]
    G = [i@b_vec for i in indices]
    incl_G = np.ones(len(G))

    # Determine which G vectors lie on BZ boundary
    #   If G/2 is closer to Gamma than to another reciprocal lattice vector
    for i,g1 in enumerate(G):
      half_G = g1/2
      for j,g2 in enumerate(G[:len(G)//2]):
        if i != j and norm(half_G) >= norm(half_G-g2):
          incl_G[i] = 0
    planes = np.array([g for i,g in enumerate(G) if incl_G[i]])

    corners = []
    sqr = lambda v : np.sum(v*v)
    # Search combinations of planes and save intersections as corners
    for i,p1 in enumerate(planes):
      for j,p2 in enumerate(planes[i+1:]):
        for k,p3 in enumerate(planes[j+1:]):
          M = np.array([p1,p2,p3])
          magG = .5 * np.array([sqr(p1),sqr(p2),sqr(p3)])
          if not np.isclose(det(M),0.):
            c = solve(M,magG)
            corners.append(c)

    corners = np.array(corners)
    # Set near zero vlues to zero explicitly
    for i,c in enumerate(corners):
      for j,v in enumerate(c):
        if np.isclose(v,0):
          corners[i][j] = 0.

    # Select corners closer to Gamma than to another reciprocal lattice vector
    # Then, eliminate any repeated corners
    incl_C = np.ones(len(corners))
    for i,c in enumerate(corners):
      for g in G:
        if norm(c) > norm(c-g)+1e-10:
          incl_C[i] = 0

    # Eliminate duplicate corners
    corners = np.unique(np.array([t for i,t in enumerate(corners) if incl_C[i]]), axis=0)

    lines = []
    # Compute lines of points which share two planes
    for ci,c1 in enumerate(corners):
      for c2 in corners[ci+1:]:
        for pi,p1 in enumerate(planes):
          for p2 in planes[pi+1:]:
            incl = True
            for p in [(p1,c1), (p1,c2), (p2,c1), (p2,c2)]:
              d = np.sum(p[0]*(p[1] - p[0]/2))
              if not np.isclose(d, 0.):
                incl = False
                break
            if incl:
              lines.append((c1,c2))

    # Adjust planes to represent midpoints between reciprocal lattice points
    planes = np.array(planes)/2

    self.frame = actor.streamtube(lines, colors=(1,1,1), linewidth=self.frame_width)
    self.scene.add(self.frame)

