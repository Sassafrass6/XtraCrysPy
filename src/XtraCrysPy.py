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
    print(x,y)
    self.frame_checkbox.position = (10, y-35)

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


  def render_iso_surface ( self, data, iso_val=0, color=(1,.3,0) ):
    from scipy.spatial import ConvexHull
    from .iso_surface import iso_surface

    dmin = np.min(data)
    data -= dmin
    iso_val -= dmin

    scale = 255/np.max(data)
    data *= scale
    iso_val *= scale

    pts = self.bound_points
    origin = np.array([-.5,-.5,-.5])
    hull = (pts, ConvexHull(pts).simplices)

    self.surface = iso_surface(data, iso_val, origin, color, hull)
    self.scene.add(self.surface)


  def start_crystal_view ( self ):
    from fury import pick

    self.picker = pick.PickingManager()
    self.smanager.iren.AddObserver('WindowResizeEvent', self.update_buttons)
    if self.axes:
      self.smanager.iren.AddObserver('InteractionEvent', self.update_axes)
    self.smanager.start()
