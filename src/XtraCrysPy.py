from fury.utils import update_actor
import numpy as np

class XtraCrysPy:

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True, background=(0,0,0), perspective=False, image_prefix='XCP_Image', resolution=4 ):
    '''
    Arguments:
    '''
    from fury.interactor import CustomInteractorStyle
    from fury.data import read_viz_icons
    from fury import ui,window
    from .Line2D import Line2D

    self.wsize = size
    self.frame = None
    self.picker = None

    background = np.array(background)
    if background.max() > 1:
      background = np.array(background) / 255
    bg_norm = np.linalg.norm(background)
    self.font_color = (0,0,0) if bg_norm > 0.5 else (1,1,1)

    self.scene = window.Scene()
    self.scene.background(background)
    self.smanager = window.ShowManager(self.scene, size=size, order_transparent=True)
    self.smanager.initialize()

    self.shift_step = 0.25
    self.ui_visible = True
    if not perspective:
      self.scene.projection('parallel')

    checkbox = ['Boundary']
    initial = checkbox if boundary else []
    self.frame_checkbox = ui.Checkbox(checkbox, initial, font_size=24, font_family='Arial', position=(10,self.wsize[1]-35))
    self.scene.add(self.frame_checkbox)
    self.frame_checkbox.on_change = self.draw_frame
    for label in self.frame_checkbox.labels:
      tactor = self.frame_checkbox.options[label].text.actor
      tactor.GetTextProperty().SetColor(self.font_color)

    self.fprefix = image_prefix
    self.image_res = resolution
    cam_icon = ('camera',read_viz_icons(fname='camera.png'))
    self.cam_panel = ui.Panel2D((40,40), (5, size[1]-108), opacity=0)
    self.cam_button = ui.Button2D(icon_fnames=[cam_icon], size=(40,40))
    self.cam_button.on_left_mouse_button_clicked = self.camera_engaged
    self.cam_panel.add_element(self.cam_button, (0,0))
    self.scene.add(self.cam_panel)

    self.bound_points = None

    self.surfaces = None

    self.surface_index = 0
    self.surface_slider = None

    self.axes = axes
    self.axes_marker = None
    self.axes_visible = True

    self.smanager.add_iren_callback(self.key_press_callback, event='KeyPressEvent')


  def left_click ( self, obj, event ):
    pass


  def key_press_callback ( self, obj, event ):

    key = obj.GetKeySym().lower()
    shift = obj.GetShiftKey()
    control = obj.GetControlKey()

    if control:
      if key in ['q', 'w', 'z']:
        exit()

    if shift:
      if key == 's':
        self.save_image(self.fprefix)
      elif key == 'c':
        self.camera_default_position()
        self.smanager.render()

    if key == 'u':
      self.toggle_ui()
    elif key == 'a':
      self.toggle_axes()
    elif key == 'b':
      self.toggle_frame()
      self.smanager.render()

    elif key in ['up', 'down', 'left', 'right']:

      camera = self.scene.GetActiveCamera()
      if shift:
        trans = 0
        step = self.shift_step * (1 if control else 5)
        camera.OrthogonalizeViewUp()
        up = np.array(camera.GetViewUp())
        focal = np.array(camera.GetFocalPoint())
        position = np.array(camera.GetPosition())
        if key == 'up':
          trans = -step * up
        elif key == 'down':
          trans = step * up
        elif key == 'left':
          right = np.cross(focal - position, up)
          trans = step * right / np.linalg.norm(right)
        elif key == 'right':
          right = np.cross(focal - position, up)
          trans = -step * right / np.linalg.norm(right)
        camera.SetPosition(position + trans)
        camera.SetFocalPoint(focal + trans)

      else:
        step = 1 if control else 5
        if key == 'up':
          camera.Elevation(-step)
        elif key == 'down':
          camera.Elevation(step)
        elif key == 'left':
          camera.Azimuth(step)
        elif key == 'right':
          camera.Azimuth(-step)

      camera.OrthogonalizeViewUp()
      self.scene.ResetCameraClippingRange()
      self.smanager.render()


  def camera_default_position ( self ):
    cam = self.scene.GetActiveCamera()
    cam.SetPosition(self.cam_defaults[0])
    cam.SetFocalPoint(self.cam_defaults[1])
    cam.SetViewUp(self.cam_defaults[2])
    self.scene.ResetCamera()


  def camera_engaged ( self, iren, caller, event ):
    self.save_image(self.fprefix)


  def save_image ( self, fprefix ):
    from fury.lib import RenderWindowInteractor,WindowToImageFilter,numpy_support
    from fury.io import save_image
    from os.path import isfile

    print('Taking snapshot...')
    window_to_image_filter = WindowToImageFilter()
    window_to_image_filter.SetInput(self.scene.GetRenderWindow())
    window_to_image_filter.SetScale(self.image_res)
    window_to_image_filter.Update()

    vtk_image = window_to_image_filter.GetOutput()
    h, w, _ = vtk_image.GetDimensions()
    vtk_array = vtk_image.GetPointData().GetScalars()
    components = vtk_array.GetNumberOfComponents()
    arr = numpy_support.vtk_to_numpy(vtk_array).reshape(w, h, components)

    nim = 0
    fname = lambda n : '{}.{}.png'.format(fprefix, n)
    while isfile(fname(nim)):
      nim += 1

    fname = fname(nim)
    print('  Saving image as: {}'.format(fname))
    save_image(arr, fname, compression_quality=100)

    self.smanager.render()


  def toggle_ui ( self ):

    vis = self.ui_visible = not self.ui_visible
    self.cam_panel.set_visibility(self.ui_visible)
    self.frame_checkbox.set_visibility(self.ui_visible)

    if self.surface_slider is not None and not self.disp_all_iso:
      self.surface_slider.set_visibility(vis)
      self.surface_slider.update()

    self.smanager.render()


  def toggle_axes ( self ):
    if self.axes:
      vis = self.axes_visible = not self.axes_visible
      self.axes_marker.SetEnabled(vis)
      self.smanager.render()


  def toggle_frame ( self ):
    option = self.frame_checkbox.options['Boundary']
    option.checked = not option.checked
    if option.checked:
      option.select()
    else:
      option.deselect()
    self.frame_checkbox._handle_option_change(option)
    if self.frame is not None:
      self.draw_frame(self.frame_checkbox)
      self.smanager.render()


  def draw_frame ( self, checkboxes ):
    if self.frame is not None:
      if 'Boundary' in checkboxes.checked_labels:
        self.scene.add(self.frame)
      else:
        self.scene.rm(self.frame)


  def iso_slider_handle_color (self, i_ren, _obj, _slider):
    self.surface_slider.handle.color = self.font_color
    i_ren.force_render()


  def update_buttons ( self, caller, event ):
    x,y = self.scene.GetSize()
    self.cam_panel.position = (5, y-108)
    self.frame_checkbox.position = (10, y-35)
    if self.surface_slider is not None:
      self.surface_slider.center = (x/2, y-50)


  def update_iso_surface ( self, slider ):
    sind = int(np.round(slider.value)) - 1
    if sind != self.surface_index:
      su,ar = self.surfaces[self.surface_index]
      self.scene.rm(su)
      self.scene.add(self.surfaces[sind][0])
      if ar is not None:
        self.scene.rm(ar)
      if self.surfaces[sind][1] is not None:
        self.scene.add(self.surfaces[sind][1])
      self.surface_index = sind


  def render_iso_surface ( self, lattice, origin, data, arrows, iso_vals=0, colors=(255,110,0,255), arrow_colors=(255,100,0,255), arrow_scale=0.25, arrow_anchor='mid', arrow_spacing=0.01, disp_all=False, clip_planes=None, clip_boundary=True, nsc=(1,1,1) ):
    '''
      Draw an isosurface from volumetric data. Data may be colored with the colors argument, either as a single color or with a color for each voxel. Arrows can be displayed by providing arrows with one normal for each data point. The arrows can be independently colored with arrow_colors. Additionally, the data can be clipped by specifying plane points and normals in the clip_planes argument. clip_planes must be of dimension (2,N,3) where N is an arbitrary number of planes to clip on. The first dimension specifies points on index 0 and normals on index 1.
      Arguments:
        self (XtraCrysPy):
        lattice (ndarray or list): Array similar to clip_planes (2,N,3), bounding the data
        origin (ndarray or list): Arraw specifying the translation required to bring the data's center to the origin. Dimension (3,)
        data (ndarray or list): Array of scalars on which to compute the isosurface. Dimension (X,Y,Z) for arbitrary X,Y, and Z
        arrows (ndarray or list): Array of normals for arrows to draw on the surface. Dimension (X,Y,Z,3) with X,Y,Z determined by data dimensions
        iso_vals (float or list): Value or list of iso-values to generate surfaces
        colors (ndarray or list): Colors for the surface. Can specify for each isovalue and for each voxel, or just choose single colors. Accepted dimensions are (A,), (N,A), (X,Y,Z,A), or (N,X,Y,Z,A), where X,Y,Z are determined by data dimensions, N is the number of isovalues provided, and A is either 3 or 4 for RGB or RGBA colors respectively
        arrow_colors (ndarray or list): Colors for the arrows, same specifications as the surface colors
        arrow_anchor (str): Anchor position for arrows. Options are 'mid', 'tip', and 'tail'
        arrow_spacing (float): Tolerance for cleaning how many arrows can appear in a certain region. Increase this value to recuce the arrow density.
        disp_all (bool): True draws all surfaces at once, False adds a slider for choosing displayed surface.
        clip_planes (ndarray or list): Specify plane points and normals for cutting the isosurface and arrows. Dimension (2,N,3) where N is an arbitrary number of planes to clip on. The first dimension specifies points on index 0 and normals on index 1.
        clip_boundary (bool): Setting False disables clipping of the isosurface within the first BZ.
    '''
    from scipy.spatial import ConvexHull
    from .iso_surface import iso_surface
    from fury import ui

    self.disp_all_iso = disp_all
    data = np.array(data)
    if len(data.shape) != 3:
      raise Exception('Argument data must be a 3D array')
    else:
      data = data[::-1, ::-1, ::-1]

    dmin = np.min(data)
    data -= dmin
    iso_vals -= dmin

    dmax = np.max(data)
    scale = 255/dmax
    data *= scale
    iso_vals *= scale

    print('Data range: [{}, {}]'.format(dmin, dmax+dmin))

    if isinstance(iso_vals, (float,int)):
      iso_vals = np.array([iso_vals], dtype=float)
    else:
      iso_vals = np.array(iso_vals, dtype=float)

    if len(iso_vals.shape) != 1:
      raise Exception('Argumenta iso_vals must be a number or a 1D array of numbers')

    niv = iso_vals.shape[0]
    colors = np.array(colors, dtype=float)
    arrow_colors = np.array(arrow_colors, dtype=float)

    def format_colors_array ( col ):
      cshape = col.shape
      one_col = True
      cdim = len(cshape)

      if np.max(col) <= 1:
        col *= 255

      if cdim == 1:
        if cshape[0] not in [3, 4]:
          raise Exception('Colors should be a tuple (R,G,B) or (R,G,B,A) with values in the range [0,255]')
        if cshape[0] == 3:
          col = np.concatenate([col,[255]])
        col = np.array([col] * niv)
      elif cdim == 2:
        if cshape[0] != niv or cshape[1] not in [3, 4]:
          raise Exception('Must provide one 3-color for each iso_val')
        if cshape[1] == 3:
          tcol = np.empty((niv,4), dtype=float)
          tcol[:,:3] = col[:]
          tcol[:,3] = 255
          col = tcol
      elif cdim == 3:
        raise Exception('Must provide colors in the same shape as data')
      elif cdim == 4 or cdim == 5:
        ci = cdim - 4
        one_col = False
        for i in range(3):
          if cshape[ci+i] != data.shape[i]:
            raise Exception('Must provide one 3-color for each data point')
        if cdim == 4:
          if cshape[-1] == 3:
            tcol = np.empty(list(cshape[:-1])+[4], dtype=float)
            tcol[:,:,:,:3] = col
            tcol[:,:,:,3] = 255
            col = tcol
          col = np.array([col[::-1,::-1,::-1]] * niv)
        else:
          if cshape[0] != niv:
            raise Exception('First array dimension must match the number of iso_vals')
          elif cshape[-1] != 3:
            raise Exception('Colors must be provided as 3-colors')
          else:
            if cshape[-1] == 3:
              tcol = np.empty(list(cshape[:-1])+[4], dtype=float)
              tcol[:,:,:,:,:3] = col
              tcol[:,:,:,:,3] = 255
              col = tcol
            col = col[:, ::-1, ::-1, ::-1, :]
      else:
        raise ValueError('Invalid dimensions for colors array')

      return col, one_col

    colors,one_col = format_colors_array(colors)
    arrow_colors,one_acol = format_colors_array(arrow_colors)

    if clip_planes is None:
      bound_planes = self.bound_planes if clip_boundary else None
    else:
      clip_planes = np.array(clip_planes)
      if clip_planes.shape[0] != 2 or clip_planes.shape[2] != 3:
        raise Exception('clip_planes must have shape (2,N,3) where N is the number of planes to clip on. The first dimension has index 0 for points and index 1 for normals.')
      if clip_boundary:
        nbp = self.bound_planes.shape[1]
        ncp = clip_planes.shape[1]
        bound_planes = np.empty((2,nbp+ncp,3))
        bound_planes[:,:nbp,:] = self.bound_planes[:,:,:]
        bound_planes[:,nbp:,:] = clip_planes[:,:,:]
      else:
        bound_planes = clip_planes

    ds = data.shape
    scshape = [(nsc[i]+1)*s for i,s in enumerate(ds)]
    ndata = np.empty(scshape, dtype=float)
    if not one_col:
      ncols = np.empty([niv]+scshape+[4], dtype=float)
    if arrows is not None:
      narr = np.empty(scshape+[3], dtype=float)
      if not one_acol:
        nacols = np.empty([niv]+scshape+[4], dtype=float)
    for i in range(nsc[0]+1):
      for j in range(nsc[1]+1):
        for k in range(nsc[2]+1):
          dw = (i, j, k)
          up = (i+1, j+1, k+1)
          dw0,dw1,dw2 = (dw[i]*ds[i] for i in range(3))
          up0,up1,up2 = (up[i]*ds[i] for i in range(3))
          ndata[dw0:up0, dw1:up1, dw2:up2] = data[:, :, :]
          if not one_col:
            for ic in range(niv):
              ncols[ic, dw0:up0, dw1:up1, dw2:up2, :] = colors[ic, :, :, :, :]
          if arrows is not None:
            narr[dw0:up0, dw1:up1, dw2:up2, :] = arrows[:, :, :, :]
            if not one_acol:
              nacols[ic, dw0:up0, dw1:up1, dw2:up2, :] = arrow_colors[ic, :, :, :, :]

    data = ndata
    ds = data.shape
    if not one_col:
      colors = ncols
    if arrows is not None:
      arrows = narr
      if not one_acol:
        arrow_colors = nacols
    for i in range(3):
      data = np.roll(data, ds[i]//4, axis=i)
      if not one_col:
        for ic in range(niv):
          colors[ic] = np.roll(colors[ic], ds[i]//4, axis=i)
      if arrows is not None:
        arrows = np.roll(arrows, ds[i]//4, axis=i)
        if not one_acol:
          arrow_colors = np.roll(arrow_colors, ds[i]//4, axis=i)

    nsurf = 0
    if self.surfaces is None:
      self.surfaces = []
    else:
      nsurf = len(self.surfaces)
    grid_spacing = [(nsc[i]+1)/s for i,s in enumerate(data.shape)]
    for i,iv in enumerate(iso_vals):
      iso = iso_surface(data, grid_spacing, iv, origin, colors[i],
                        bound_planes=bound_planes,
                        skew=lattice, arrows=arrows,
                        arrow_colors=arrow_colors[i],
                        arrow_scale=arrow_scale,
                        arrow_anchor=arrow_anchor,
                        arrow_spacing=arrow_spacing)
      self.surfaces.append(iso)

    if len(self.surfaces) > nsurf:
      if not disp_all:
        for su,ar in self.surfaces[:nsurf]:
          self.scene.rm(su)
          if ar is not None:
            self.scene.rm(ar)
        if len(self.surfaces) == 1:
          self.scene.add(self.surfaces[0][0])
          if self.surfaces[0][1] is not None:
            self.scene.add(self.surfaces[0][1])
        else:
          if self.surface_slider is None:
            self.surface_index = 0
            self.surface_slider = ui.LineSlider2D(center=(self.wsize[0]/2,self.wsize[1]-50), initial_value=1, orientation='horizontal', min_value=1, max_value=len(self.surfaces), text_template='{value:0.0f}')
            self.surface_slider.on_change = self.update_iso_surface
            self.surface_slider.text.color = self.font_color
            self.surface_slider.handle.color = self.font_color
            self.surface_slider.handle.on_left_mouse_button_released = self.iso_slider_handle_color
            self.scene.add(self.surface_slider)
          else:
            self.surface_slider.max_value = len(self.surfaces)
            self.surface_slider.set_visibility(True)
            self.surface_slider.update()
          self.scene.add(self.surfaces[self.surface_index][0])
          if self.surfaces[self.surface_index][1] is not None:
            self.scene.add(self.surfaces[self.surface_index][1])
      else:
        if self.surface_slider is not None:
          self.surface_slider.set_visibility(False)
        for su,ar in self.surfaces:
          self.scene.add(su)
          if ar is not None:
            self.scene.add(ar)


  def start_crystal_view ( self ):
    '''
      Begin the render sequence and allow interaction
    '''
    from fury import pick

    cam = self.scene.GetActiveCamera()
    self.cam_defaults = (cam.GetPosition(),
                         cam.GetFocalPoint(),
                         cam.GetViewUp())

    self.picker = pick.PickingManager()
    self.smanager.iren.AddObserver('WindowResizeEvent', self.update_buttons)

    if self.axes:
      from vtk import vtkAxesActor,vtkOrientationMarkerWidget
      self.axes_marker = vtkOrientationMarkerWidget()
      self.axes_marker.SetInteractor(self.scene.GetRenderWindow().GetInteractor())
      self.axes_marker.SetOrientationMarker(vtkAxesActor())
      self.axes_marker.SetViewport(0,0,.2,.2)
      self.axes_marker.SetEnabled(1)

    self.camera_default_position()
    self.smanager.render()
    self.smanager.start()
