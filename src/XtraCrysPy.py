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


  def render_iso_surface ( self, data, iso_vals=0, colors=(1,.3,0) ):
    from scipy.spatial import ConvexHull
    from .iso_surface import iso_surface
    from fury import ui

    if isinstance(iso_vals, (float,int)):
      iso_vals = np.array([iso_vals])

    colors = np.array(colors)
    if len(colors.shape) == 1:
      colors = [colors] * len(iso_vals)
    elif len(colors.shape) == 2:
      if colors.shape[0] != len(iso_vals) or colors.shape[1] != 3:
        print('Must provide one 3-color for each iso_val')
        return
    elif len(colors.shape) == 4:
      colors = [colors] * len(iso_vals)
    elif len(colors.shape) == 5:
      print('Colors 5 shape not implemented')
      return
    else:
      print('Invalid dimensions for colors array')
      return

    dmin = np.min(data)
    data -= dmin
    iso_vals -= dmin

    scale = 255/np.max(data)
    data *= scale
    iso_vals *= scale

    hull = None
    if self.bound_points is not None:
      pts = self.bound_points
      hull = (pts, ConvexHull(pts).simplices)

    #s = data.shape
    #orig = np.array([40,40,30])
    #color = np.zeros((s[0],s[1],s[2],3), dtype=float, order='C')
    #for i,x in enumerate(color):
    #  for j,y in enumerate(x):
    #    for k,z in enumerate(y):
    #      #d = np.sqrt(np.sum((np.array([i,j,k])-orig)**2))
    #      color[i,j,k,:] = np.array([255*i/s[0],0,100]).astype(int)
          #color[i,#j,k,:] = (255*np.round(np.array([i/s[0],j/s[1],k/s[2]]),2)).astype(int)
    #color[:40,:,:,0] = 255
    #color[40:,:,:,2] = 255
    #color[:20,:,:,0] = 255
    #color[:,:20,:,1] = 255
    #color[:,:,:20,2] = 255
    #color[:40,:,:,0] = 255
    #color[40:,:,:,2] = 255
    #color = np.roll(color, 40, axis=2)
    #color = np.roll(color, 20, axis=0)
    #color[30:50,:,:,:] = color[49:29:-1,:,:,:]
    #color[:,:,:,:] = color[::-1,:,:,:]
    #color = np.roll(color, -20, axis=1)
    #color = np.roll(color, -20, axis=2)
    #color[:,:40,:,:] = cc[:,:39:-1,:,:]
    #color[:,40:,:,:] = cc[:,:40,:,:]
    #color = np.swapaxes(color, 0, 2)
    #color = np.swapaxes(np.swapaxes(color, 0, 1), 0, 2)

    self.surfaces = []
    origin = np.array([-.5,-.5,-.5])
    for i,iv in enumerate(iso_vals):
      self.surfaces.append(iso_surface(data, iv, origin, colors[i], hull))
    self.scene.add(self.surfaces[0])

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
