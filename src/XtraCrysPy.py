from fury.utils import update_actor
import numpy as np

class XtraCrysPy:

  def __init__ ( self, size=(1024, 1024), axes=True, perspective=False ):
    '''
    Arguments:
    '''
    from .Line2D import Line2D
    from fury import ui,window

    self.wsize = size
    self.frame = None
    self.picker = None

    self.scene = window.Scene()
    #self.scene.projection('parallel')
    self.smanager = window.ShowManager(self.scene, size=size, order_transparent=True)
    self.smanager.initialize()

    checkbox = ['Boundary']
    self.frame_checkbox = ui.Checkbox(checkbox, checkbox, font_size=24, font_family='Arial', position=(10,self.wsize[1]-35))
    self.scene.add(self.frame_checkbox)
    self.frame_checkbox.on_change = self.toggle_frame

    self.bound_points = None

    self.surface_index = 0
    self.surface_slider = None

    self.axes = axes
    if axes:
      try:
        self.axis_eles = [None]*3
        self.axis_lines = [None]*3
        self.axis_width = ax_wid = 200
        self.ax_panel = ui.Panel2D(size=(ax_wid,ax_wid), color=(0,0,0))
        self.ax_panel.center = (ax_wid/2, ax_wid/2)
        self.scene.add(self.ax_panel)

        axis_vecs = np.array([[1,0,0],[0,1,0],[0,0,1]])
        centers = self.axis_endpoint_positions(axis_vecs)
        for i in range(3):
          self.axis_eles[i] = ui.Disk2D(outer_radius=8, center=centers[i], color=axis_vecs[i])
          self.axis_lines[i] = Line2D((ax_wid/2,ax_wid/2), centers[i], color=axis_vecs[i])
          self.axis_lines[i].width = 5

          self.scene.add(self.axis_eles[i])
          self.scene.add(self.axis_lines[i])
      except Exception as e:
        print(e)
        print('Could not import Line2D from fury. Instal Line2D branch from Sassafrass6 GitHub')
        print('Will not display coordinate axes')
        self.axes = False


  def axis_endpoint_positions ( self, vecs ):
    import numpy as np
    centers = np.empty((3,3), dtype=float)
    for i in range(3):
      wid = self.axis_width / 2
      centers[i,:] = np.array([wid,wid,0]) + 80*vecs[i]
    return centers[:,:2]


  def left_click ( self, obj, event ):
    pass


  def toggle_frame ( self, checkboxes ):
    if self.frame is not None:
      if 'Boundary' in checkboxes.checked_labels:
        self.scene.add(self.frame)
      else:
        self.scene.rm(self.frame)


  def update_buttons ( self, caller, event ):
    x,y = self.scene.GetSize()
    self.frame_checkbox.position = (10, y-35)
    if self.surface_slider is not None:
      self.surface_slider.center = (x/2, y-50)


  def update_axes ( self, caller, event ):

    camera = self.scene.GetActiveCamera()

    y = np.array(camera.GetViewUp())
    y /= np.linalg.norm(y)
    z = np.array(camera.GetFocalPoint()) - np.array(camera.GetPosition())
    z /= np.linalg.norm(z)
    x = np.cross(y, z)
    x /= np.linalg.norm(x)

    rot = []
    unit = np.eye(3)
    cam_vecs = np.array([-x,y,z])
    for i,u in enumerate(unit[:2]):
      v = np.cross(cam_vecs[i], u)
      r = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
      rot.append(np.eye(3) + r + r@r/(1+cam_vecs[i].dot(u)))
      cam_vecs = rot[-1] @ cam_vecs.T

    axis_vecs = np.eye(3) @ np.linalg.inv(rot[0]) @ np.linalg.inv(rot[1])

    for i,v in enumerate(axis_vecs):
        axis_vecs[i] /= np.linalg.norm(v)
    centers = self.axis_endpoint_positions(axis_vecs)
    axis_order = np.argsort([v[2] for v in axis_vecs])
    for i,o in enumerate(axis_order):
      color = np.zeros(3); color[o] = 1
      self.axis_lines[i].color = color
      self.axis_lines[i].p2 = centers[o]
      self.axis_eles[i].color = color
      self.axis_eles[i].center = centers[o]


  def update_iso_surface ( self, slider ):
    sind = int(np.round(slider.value)) - 1
    if sind != self.surface_index:
      self.scene.rm(self.surfaces[self.surface_index])
      self.scene.add(self.surfaces[sind])
      self.surface_index = sind


  def render_iso_surface ( self, data, iso_vals=0, colors=(255,110,0) ):
    from scipy.spatial import ConvexHull
    from .iso_surface import iso_surface
    from fury import ui

    data = np.array(data)
    if len(data.shape) != 3:
      raise Exception('Argument data must be a 3D array')

    if isinstance(iso_vals, (float,int)):
      iso_vals = np.array([iso_vals], dtype=float)
    else:
      iso_vals = np.array(iso_vals, dtype=float)

    if len(iso_vals.shape) != 1:
      raise Exception('Argumenta iso_vals must be a number or a 1D array of numbers')

    colors = np.array(colors, dtype=float)
    cshape = colors.shape

    cdim = len(cshape)
    niv = iso_vals.shape[0]
    colors = np.array(colors)
    if cdim == 1:
      if cshape[0] != 3:
        raise Exception('Colors should be a tuple (R,G,B) with values in the range [0,255]')
      colors = np.array([colors] * niv)
    elif cdim == 2:
      if cshape[0] != niv or cshape[1] != 3:
        raise Exception('Must provide one 3-color for each iso_val')
    elif cdim == 4 or cdim == 5:
      ci = cdim - 4
      for i in range(3):
        if cshape[ci+i] != data.shape[i]:
          raise Exception('Must provide one 3-color for each data point')
      if cdim == 4:
        colors = np.array([colors] * niv)
    elif cdim == 5 and cshape[0] != niv:
      raise Exception('First array dimension must match the number of iso_vals')
    else:
      raise ValueError('Invalid dimensions for colors array')

    if np.max(colors[0]) <= 1:
      for i in range(niv):
        colors[i] *= 255

    dmin = np.min(data)
    data -= dmin
    iso_vals -= dmin

    scale = 255/np.max(data)
    data *= scale
    iso_vals *= scale

    if len(colors.shape)  > 2:

      resize = False
      cshape = colors.shape
      mdim = np.max(data.shape)
      for i in range(1,4):
        if cshape[i] != mdim:
          resize = True

      if resize:
        cs = cshape[1:]
        nshape = [niv, mdim, mdim, mdim, 3]
        ncol = np.zeros(nshape, dtype='uint8')
        print('Padding colors with zeros to size {}x{}x{}'.format(*nshape[1:]))
        for i in range(niv):
          ncol[i,:cs[0],:cs[1],:cs[2],:] = colors[i][:,:,:,:]
        cs = nshape
        colors = ncol
        ncol = np.zeros([niv]+[2*s for s in nshape[1:-1]]+[3])
        for i in range(2):
          for j in range(2):
            for k in range(2):
              dw = (i, j, k)
              up = (i+1, j+1, k+1)
              for iv in range(niv):
                ncol[iv, dw[0]*cs[1]:up[0]*cs[1], dw[1]*cs[2]:up[1]*cs[2], dw[2]*cs[3]:up[2]*cs[3]] = colors[iv,:,:,:,:]
        del colors
        colors = ncol
        for i in range(3):
          colors = np.roll(colors, data.shape[i]//2, axis=i+1)

    hull = None
    if self.bound_points is not None:
      pts = self.bound_points
      hull = (pts, ConvexHull(pts).simplices)

    ds = data.shape
    ndata = np.empty([2*s for s in ds])
    for i in range(2):
      for j in range(2):
        for k in range(2):
          dw = (i, j, k)
          up = (i+1, j+1, k+1)
          ndata[dw[0]*ds[0]:up[0]*ds[0], dw[1]*ds[1]:up[1]*ds[1], dw[2]*ds[2]:up[2]*ds[2]] = data[:,:,:]

    data = ndata
    for i in range(3):
      data = np.roll(data, data.shape[i]//4, axis=i)

    self.surfaces = []
    origin = 2*np.array([-.5,-.5,-.5])
    for i,iv in enumerate(iso_vals):
      self.surfaces.append(iso_surface(data, iv, origin, colors[i], bound_polys=hull, lattice=self.rlattice))
    self.scene.add(self.surfaces[0])

    if len(self.surfaces) > 1:
      self.surface_index = 0
      self.surface_slider = ui.LineSlider2D(center=(self.wsize[0]/2,self.wsize[1]-50), initial_value=1, orientation='horizontal', min_value=1, max_value=iso_vals.shape[0])
      self.surface_slider.on_change = self.update_iso_surface
      self.scene.add(self.surface_slider)


  def start_crystal_view ( self ):
    from fury import pick

    self.picker = pick.PickingManager()
    self.smanager.iren.AddObserver('WindowResizeEvent', self.update_buttons)
    if self.axes:
      self.smanager.iren.AddObserver('InteractionEvent', self.update_axes)
    self.smanager.start()
