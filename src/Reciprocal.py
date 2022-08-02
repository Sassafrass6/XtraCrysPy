from .XtraCrysPy import XtraCrysPy
from .Model import Model
from fury import actor
import numpy as np

class Reciprocal ( XtraCrysPy ):

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True, background=(0,0,0), perspective=False, model=None, frame_width=4e-3, image_prefix='XCP_Image', resolution=4 ):
    super().__init__(size, axes, boundary, background,
                     perspective, image_prefix, resolution)

    self.model = model
    self.point_actors = []
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

      elif isinstance(model, str):
        from os.path import isfile
        if not isfile(model):
          raise FileNotFoundError('File {} not found'.format(model))
        self.model = Model(fname=model)
        self.rlattice = self.model.rlattice

      elif isinstance(model, Model):
        self.rlattice = model.rlattice

      else:
        raise ValueError('Argument \'model\' must be a Model object or a dictionary.')

      self.show_k_axes = True
      zeros = np.zeros((3,3), dtype=float)
      self.k_axes = actor.arrow(zeros, self.rlattice, zeros,
                                tip_radius=.005, tip_length=0.025,
                                shaft_radius=0.002, repeat_primitive=False)
      self.scene.add(self.k_axes)

      self.shift_step = 0.01
      self.bravais_boundaries(render=boundary)
      cam = self.scene.GetActiveCamera()
      self.smanager.render()
      self.cam_defaults = (cam.GetPosition(),
                           cam.GetFocalPoint(),
                           cam.GetViewUp())


  def toggle_k_axes ( self ):
    self.show_k_axes = not self.show_k_axes
    if self.k_axes is not None:
      if self.show_k_axes:
        self.scene.add(self.k_axes)
      else:
        self.scene.rm(self.k_axes)
      self.smanager.render()


  def key_press_callback ( self, obj, event ):

    key = obj.GetKeySym().lower()
    shift = obj.GetShiftKey()
    control = obj.GetControlKey()

    if key == 'k':
      self.toggle_k_axes()

    else:
      super().key_press_callback(obj, event)


  def display_points ( self, points, colors=(1,1,1), radii=0.04 ):
    from fury import actor

    points = np.array(points)

    colors = np.array(colors)
    if len(colors.shape) == 1:
      colors = np.array([colors]*points.shape[0])
    if colors.shape != (points.shape[0], 3):
      raise ValueError('Argument colors should be either single valued or the same shape as points.')

    if isinstance(radii, (float,int)):
      radii = np.array([radii]*points.shape[0])
    else:
      radii = np.array(radii)
      if len(radii.shape) != 1 or radii.shape[0] != points.shape[0]:
        raise ValueError('Argument radii should be either single valued or have one radius for each point.')

    try:
      p_actors = actor.sphere(points, colors, radii, use_primitive=False)
    except:
      p_actors = actor.sphere(points, colors, radii)

    self.scene.add(p_actors)
    self.point_actors.append(p_actors)


  def bravais_boundaries ( self, render=True ):
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
    self.bound_planes = np.array([np.array(planes)/2]*2)
    self.bound_points = np.unique(np.array(lines).reshape((2*len(lines),3)), axis=0)

    self.frame = actor.streamtube(lines, colors=(1,1,1), linewidth=self.frame_width)
    if render:
      self.scene.add(self.frame)

    self.scene.ResetCamera()


  def render_iso_surface ( self, data, origin=(0,0,0), arrows=None, iso_vals=0, colors=(255,110,0), arrow_colors=(255,100,0), arrow_scale=0.25, arrow_anchor='mid', arrow_spacing=0.01, disp_all=False, clip_planes=None, clip_boundary=True ):
    '''
      Draw an isosurface from volumetric data. Data may be colored with the colors argument, either as a single color or with a color for each voxel. Arrows can be displayed by providing arrows with one normal for each data point. The arrows can be independently colored with arrow_colors. Additionally, the data can be clipped by specifying plane points and normals in the clip_planes argument. clip_planes must be of dimension (2,N,3) where N is an arbitrary number of planes to clip on. The first dimension specifies points on index 0 and normals on index 1.
      Arguments:
        self (XtraCrysPy):
        data (ndarray or list): Array of scalars on which to compute the isosurface. Dimension (X,Y,Z) for arbitrary X,Y, and Z
        origin (tuple, list, or ndarray): 3-vector origin to offset the surface position manually.
        arrows (ndarray or list): Array of normals for arrows to draw on the surface. Dimension (X,Y,Z,3) with X,Y,Z determined by data dimensions
        iso_vals (float or list): Value or list of iso-values to generate surfaces
        colors (ndarray or list): Colors for the surface. Can specify for each isovalue and for each voxel, or just choose single colors. Accepted dimensions are (A,), (N,A), (X,Y,Z,A), or (N,X,Y,Z,A), where X,Y,Z are determined by data dimensions, N is the number of isovalues provided, and A is either 3 or 4 for RGB or RGBA colors respectively
        arrow_colors (ndarray or list): Colors for the arrows, same specifications as the surface colors
        arrow_scale (float): Scale for the displayed arrows
        arrow_anchor (str): Anchor position for arrows. Options are 'mid', 'tip', and 'tail'
        arrow_spacing (float): Tolerance for cleaning how many arrows can appear in a certain region. Increase this value to recuce the arrow density.
        disp_all (bool): True draws all surfaces at once, False adds a slider for choosing displayed surface.
        clip_planes (ndarray or list): Specify plane points and normals for cutting the isosurface and arrows. Dimension (2,N,3) where N is an arbitrary number of planes to clip on. The first dimension specifies points on index 0 and normals on index 1.
        clip_boundary (bool): Setting True disables clipping of the isosurface within the first BZ.
    '''
    origin = np.array(origin, dtype=float)
    origin -= 1 - 1/np.array(data.shape)
    super().render_iso_surface(self.rlattice, origin, data, arrows, iso_vals, colors, arrow_colors, arrow_scale, arrow_anchor, arrow_spacing, disp_all, clip_planes, clip_boundary)

