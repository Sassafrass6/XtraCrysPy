from fury.utils import colors_from_actor, update_actor
from fury.utils import vertices_from_actor
from .XtraCrysPy import XtraCrysPy
import numpy as np

class Atomic ( XtraCrysPy ):

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True,
                 background=(0,0,0), perspective=False, model=None,
                 params={}, multi_frame=False, nsc=(1,1,1),
                 bond_type='Stick', sel_type='Chain', unit='angstrom',
		 runit='degree', constrain_atoms=False,
                 image_prefix='XCP_Image', resolution=4 ):
    super().__init__(size, axes, boundary, background,
                     perspective, image_prefix, resolution)
    from .Model import Model

    self.nsc = list(nsc)
    self.sel_forward = True
    self.sel_type = sel_type
    self.bond_type = bond_type
    self.constrain_atoms = constrain_atoms

    self.sel_inds = []
    self.sel_cols = []
    self.sel_bnds = []
    self.scolor = np.array((0,210,210))

    if model is not None:
      if isinstance(model, Model):
        pass
      elif isinstance(model, str):
        from os.path import isfile
        if not isfile(model):
          raise FileNotFoundError('File {} not found'.format(model))
        model = Model(params, fname=model, multi_frame=multi_frame)
      elif isinstance(model, (list, tuple)):
        from os.path import isfile
        for mfn in model:
          if not isfile(mfn):
            raise FileNotFoundError('File {} not found'.format(mfn))
        model = Model(params, fname=model, multi_frame=True)
      else:
        s = 'Argument \'model\' must be a file name, a list of file names, or a Model object.'
        raise TypeError(s)
    elif params:
      if 'lattice' not in params:
        if self.boundary:
          self.toggle_frame()
      model = Model(params, fname=None, multi_frame=multi_frame)
    else:
      model = Model(params)

    self.atoms = None
    self.frame = None
    self.model = model
    self.frame_index = 0
    self.relax = model.relax
    self.nrelax = self.model.atoms.shape[0]
    self.relax_boundary = self.model.lattice.shape[0] > 1

    self.shift_step = 0.25
    self.units = unit.lower()
    self.runits = runit.lower()
    if self.units not in ['angstrom', 'bohr']:
      print('Supported length units are angstrom and bohr')
      self.units = 'angstrom'
    if self.runits not in ['degree', 'radian']:
      print('Supported radial units are degree and radian')
      self.runits = 'degree'

    self.setup_ui()

    self.render_atomic_model()
    self.scene.ResetCamera()
    self.smanager.render()


  def setup_ui ( self ):
    from fury.data import read_viz_icons
    from fury import ui

    size = self.wsize
    icons = [('left',read_viz_icons(fname='circle-left.png')), 
              ('down',read_viz_icons(fname='circle-down.png')),
              ('right',read_viz_icons(fname='circle-right.png'))]

    self.ncell_button = ui.Button2D(icon_fnames=icons[2:0:-1],
                                    size=(50,50))
    self.sel_type_button = ui.Button2D(icon_fnames=icons[2:0:-1],
                                       size=(50,50))
    self.ncell_button.on_left_mouse_button_clicked = self.toggle_ncell_menu
    self.sel_type_button.on_left_mouse_button_clicked = self.toggle_sel_menu

    sel_types = ['Chain', 'Info', 'Angle', 'Distance']
    if self.sel_type not in sel_types:
      s = '''{} is not a valid selection type. Choose from:
               Info, Angle Chain, or Distance.'''
      self.sel_type = 'Chain'
      print(s)

    self.sel_type_menu = ui.ListBox2D(sel_types, multiselection=False,
                         font_size=24, line_spacing=2, size=(150,220),
                         background_opacity=.9)
    self.sel_type_menu.panel.color = (.1,.1,.1)

    sel = sel_types.index(self.sel_type)
    self.sel_type_menu.select(self.sel_type_menu.slots[sel])
    self.sel_type_menu.on_change = self.update_selection_type

    self.sel_menu_vis = False
    self.sel_type_menu.set_visibility(False)

    self.ncell_panel_vis = False
    self.ncell_panel = ui.Panel2D(size=(180,140), color=(.1,.1,.1), opacity=.9)
    self.ncell_panel.center = [260, size[1]-80]

    nsc = self.nsc
    nmax = [max(n, 4) for n in nsc]
    self.slider_nx = ui.LineSlider2D(min_value=1,
                                     max_value=nmax[0],
                                     initial_value=nsc[0],
                                     text_template="{value:0.0f}",
                                     font_size=12,
                                     length=140)
    self.slider_ny = ui.LineSlider2D(min_value=1,
                                     max_value=nmax[1],
                                     initial_value=nsc[1],
                                     text_template="{value:0.0f}",
                                     font_size=12,
                                     length=140)
    self.slider_nz = ui.LineSlider2D(min_value=1,
                                     max_value=nmax[2],
                                     initial_value=nsc[2],
                                     text_template="{value:0.0f}",
                                     font_size=12,
                                     length=140)
    self.slider_nx.on_change = self.update_nsc_x
    self.slider_ny.on_change = self.update_nsc_y
    self.slider_nz.on_change = self.update_nsc_z
    self.ncell_panel.add_element(self.slider_nx, (.1,.75))
    self.ncell_panel.add_element(self.slider_ny, (.1,.45))
    self.ncell_panel.add_element(self.slider_nz, (.1,.15))
    self.ncell_panel.set_visibility(False)
    self.scene.add(self.ncell_panel)

    self.sel_panel = ui.Panel2D((110,50), (10, size[1]-70), opacity=0)
    self.sel_panel.add_element(self.ncell_button, (60,0))
    self.sel_panel.add_element(self.sel_type_button, (0,0))
    self.sel_panel.add_element(self.sel_type_menu, (0,-230))
    self.scene.add(self.sel_panel)

    self.sel_text_vis = False
    self.sel_tpanel = ui.Panel2D(size=(300,80), color=(0,0,0))
    self.sel_tpanel.center = (size[1]-160, 50)
    self.sel_text = ui.TextBlock2D(font_size=30, justification='left',
                    text='', vertical_justification='middle')
    self.sel_tpanel.add_element(self.sel_text, (25,25))
    self.scene.add(self.sel_tpanel)
    self.sel_tpanel.set_visibility(False)

    if self.relax:

      left_button = ui.Button2D(icon_fnames=[icons[0]], size=(50,50))
      left_button.on_left_mouse_button_clicked = self.relax_backward

      right_button = ui.Button2D(icon_fnames=[icons[2]], size=(50,50))
      right_button.on_left_mouse_button_clicked = self.relax_forward

      self.relax_panel = ui.Panel2D((110,50),
                         (size[0]-120,size[1]-70), (0,0,0), opacity=0)
      self.relax_panel.add_element(left_button, (0,0))
      self.relax_panel.add_element(right_button, (60,0))

      stext = '1/{}'.format(self.nrelax)
      self.relax_text = ui.TextBlock2D(text=stext, 
                        justification='center', font_size=30,
                        vertical_justification='top')
      self.relax_text.actor.GetTextProperty().SetColor(self.font_color)
      self.relax_panel.add_element(self.relax_text, (55,-10))
      self.scene.add(self.relax_panel)


  def key_press_callback ( self, obj, event ):

    key = obj.GetKeySym().lower()
    shift = obj.GetShiftKey()
    control = obj.GetControlKey()

    if not shift and key == 'c':
      self.toggle_constrain()
    if not shift and key == 's':
      self.toggle_sel_menu(None,None,None)
    if not shift and key == 'n':
      self.toggle_ncell_menu(None,None,None)

    elif key in ['less', 'greater', 'comma', 'period']:
      if self.relax:
        step = 1 if not control else int(np.round(self.nrelax/20))
        if key in ['less', 'comma']:
          self.relax_backward(None, obj, event, step=step)
        else:
          self.relax_forward(None, obj, event, step=step)

    else:
      super().key_press_callback(obj, event)


  def update_buttons ( self, caller, event ):
    super().update_buttons(caller, event)
    x,y = self.scene.GetSize()
    self.sel_tpanel.center = (x-160, 50)
    self.sel_panel.position = (10, y-70)
    if self.relax:
      self.relax_panel.position = (x-120, y-60)


  def update_selection_type ( self ):
    self.clear_selection_text()
    self.sel_type = self.sel_type_menu.selected[0] 
    for _ in range(len(self.sel_bnds)):
      self.pop_sbond()
    if self.atoms is not None:
      colors = colors_from_actor(self.atoms)
      nvert = int(vertices_from_actor(self.atoms).shape[0]/self.natoms)
      for i in self.sel_inds.copy():
        self.pop_atom(colors, i, nvert)
      update_actor(self.atoms)


  def clear_selection_text ( self ):
    self.sel_text_vis = False
    self.sel_text.message = ''
    self.sel_tpanel.set_visibility(False)


  def update_relax_text ( self ):
    text = '{}/{}'.format(self.frame_index+1, self.nrelax)
    self.relax_text.message = text


  def update_selection_text ( self, text ):
    self.sel_text_vis = True
    self.sel_text.message = text
    self.sel_tpanel.set_visibility(True)


  def toggle_constrain ( self ):
    self.constrain_atoms = not self.constrain_atoms
    self.redraw_atomic_model()


  def update_nsc_x ( self, slider ):
    self.update_selection_type()
    ind = int(np.round(slider.value))
    if ind != self.nsc[0]:
      self.nsc[0] = ind
      self.redraw_atomic_model()


  def update_nsc_y ( self, slider ):
    self.update_selection_type()
    ind = int(np.round(slider.value))
    if ind != self.nsc[1]:
      self.nsc[1] = ind
      self.redraw_atomic_model()


  def update_nsc_z ( self, slider ):
    self.update_selection_type()
    ind = int(np.round(slider.value))
    if ind != self.nsc[2]:
      self.nsc[2] = ind
      self.redraw_atomic_model()


  def toggle_ui ( self ):
    super().toggle_ui()

    sel_vis = self.ui_visible and self.sel_menu_vis
    self.sel_type_menu.set_visibility(sel_vis)

    ncell_vis = self.ui_visible and self.ncell_panel_vis
    self.ncell_panel.set_visibility(ncell_vis)

    sel_tvis = self.ui_visible and self.sel_text_vis
    self.sel_tpanel.set_visibility(sel_tvis)

    self.ncell_button.set_visibility(self.ui_visible)
    self.sel_type_button.set_visibility(self.ui_visible)

    if self.relax:
      self.relax_panel.set_visibility(self.ui_visible)

    self.smanager.render()


  def toggle_sel_menu ( self, iren, caller, event ):
    self.sel_menu_vis = not self.sel_menu_vis
    self.sel_type_menu.set_visibility(self.sel_menu_vis)
    self.sel_type_button.next_icon()
    self.smanager.render()


  def toggle_ncell_menu ( self, iren, caller, event ):
    self.ncell_panel_vis = not self.ncell_panel_vis
    self.ncell_panel.set_visibility(self.ncell_panel_vis)
    self.ncell_button.next_icon()
    self.smanager.render()


  def angle ( self, ai0, ai1, ai2 ):
    from numpy.linalg import norm

    v1 = self.aposs[ai0] - self.aposs[ai1]
    v2 = self.aposs[ai2] - self.aposs[ai1]
    conv = 1 if self.runits == 'radian' else 180/np.pi
    return conv * np.arccos((v1.dot(v2))/(norm(v1)*norm(v2)))


  def distance ( self, ai1, ai2 ):
    from .conversion import BOHR_ANG
    dist = np.sqrt(np.sum((self.aposs[ai2]-self.aposs[ai1])**2))
    return dist * (BOHR_ANG if self.units=='angstrom' else 1)


  def set_atom_color ( self, mem, index, nvert, col ):
    for i in range(index*nvert, (index+1)*nvert):
      mem[i] = np.array(col)


  def pop_sbond ( self, ind=-1 ):
    tbond = self.sel_bnds.pop(ind)
    self.scene.rm(tbond)


  def push_sbond ( self ):
    from .cylinder import cylinder
    from fury import actor

    if self.sel_forward:
      ai1,ai2 = self.sel_inds[-2:]
    else:
      ai1,ai2 = self.sel_inds[:2]
      
    conn = self.aposs[ai1] - self.aposs[ai2]
    dist = np.linalg.norm(conn)
    cent = (self.aposs[ai1] + self.aposs[ai2]) / 2

    brad = self.model.bond_radius(dist, ai1, ai2, self.bond_type)
    if self.bond_type != 'Sphere':
      brad *= 0.51
      tbond = cylinder([cent], [conn], [self.scolor/255],
                             radius=brad, heights=dist, resolution=20)
    else:
      ends = [[self.aposs[ai2], self.aposs[ai1]]]
      tbond = actor.streamtube(ends, colors=self.scolor,
                               linewidth=brad)

    self.scene.add(tbond)
    if self.sel_forward:
      self.sel_bnds.append(tbond)
    else:
      self.sel_bnds.insert(0, tbond)


  def pop_atom ( self, mem, index, nvert ):
    aind = self.sel_inds.index(index)
    self.sel_inds.pop(aind)
    col = self.sel_cols.pop(aind)
    self.set_atom_color(mem, index, nvert, col)


  def push_atom ( self, mem, index, nvert ):
    if index == -1:
      return False

    color = mem[index*nvert].copy()
    if self.sel_forward:
      self.sel_inds.append(index)
      self.sel_cols.append(color)
    else:
      self.sel_inds.insert(0, index)
      self.sel_cols.insert(0, color)
    self.set_atom_color(mem, index, nvert, self.scolor)
    return True


  def selection_logic ( self, colors, index, nvert ):

    if index in self.sel_inds:
      self.clear_selection_text()

      if self.sel_type == 'Info':
        self.pop_atom(colors, index, nvert)
      else:
        nsi = len(self.sel_inds)
        sind = self.sel_inds.index(index)
        npop = nsi - sind

        if nsi == 1: 
          self.pop_atom(colors, self.sel_inds[0], nvert)
        elif nsi == 2:
          self.pop_atom(colors, index, nvert)
          self.pop_sbond()
        else:
          mid = nsi - sind - 1
          if mid >= nsi-mid:
            rng = (0, sind+1)
            self.sel_forward = False
          else:
            rng = (sind, nsi)
            self.sel_forward = True
          for _ in range(rng[0], rng[1]):
            self.pop_atom(colors, self.sel_inds[rng[0]], nvert)
            if rng[0] == 0:
              self.pop_sbond(ind=rng[0])
            else:
              self.pop_sbond(ind=rng[0]-1)

    else:
      if self.sel_type == 'Info':
        if len(self.sel_inds) != 0:
          self.pop_atom(colors, self.sel_inds[-1], nvert)
        self.push_atom(colors, index, nvert)
        indmod = index % self.model.natoms
        spec = self.model.species[indmod]
        message = 'Atom {} ({}): {}'.format(indmod, index, spec)
        self.update_selection_text(message)
        print(message)

      elif self.sel_type == 'Distance':
        if len(self.sel_inds) == 0:
          self.push_atom(colors, index, nvert)
        elif len(self.sel_inds) == 1:
          if self.push_atom(colors, index, nvert):
            self.push_sbond()
            dist = self.distance(*self.sel_inds)
            message = '{:.4f} {}\n'.format(dist, self.units)
            self.update_selection_text(message)
            dtext = 'Distance between atoms {} and {}:'
            print(dtext.format(*self.sel_inds))
            print('\t{}'.format(message))

      elif self.sel_type == 'Angle':
        if len(self.sel_inds) == 0:
          self.push_atom(colors, index, nvert)
        elif len(self.sel_inds) in [1,2]:
          if self.push_atom(colors, index, nvert):
            self.push_sbond()
            if len(self.sel_inds) == 3:
              angle = self.angle(*self.sel_inds)
              message = '{:.4f} {}\n'.format(angle, self.runits)
              self.update_selection_text(message)
              dtext = 'Angle between atoms {}, {}, and {}:'
              print(dtext.format(*self.sel_inds))
              print('\t{}'.format(message))

      elif self.sel_type == 'Chain':
        selected = self.push_atom(colors, index, nvert)
        if selected and len(self.sel_inds) > 1:
          self.push_sbond()


  def relax_forward ( self, iren, obj, event, step=1 ):
    if self.frame_index == self.nrelax-1:
      return
    if self.frame_index + step >= self.nrelax:
      self.frame_index = self.nrelax-1
    else:
      self.frame_index += step
    self.update_relax_text()
    self.update_atomic_model()


  def relax_backward ( self, iren, obj, event, step=1 ):
    if self.frame_index == 0:
      return
    if self.frame_index - step < 0:
      self.frame_index = 0
    else:
      self.frame_index -= step
    self.update_relax_text()
    self.update_atomic_model()


  def pick_atom ( self, obj, event ):

    pos = self.picker.event_position(self.smanager.iren)
    pinfo = self.picker.pick(pos, self.smanager.scene)

    vertices = vertices_from_actor(obj)
    colors = colors_from_actor(obj, 'colors')

    nvert = int(vertices.shape[0]/self.natoms)
    index = int(np.floor(pinfo['vertex']/nvert))

    self.selection_logic(colors, index, nvert)

    update_actor(obj)


  def update_boundary ( self, planes, lines ):
    from fury.actor import streamtube
    self.bound_planes = planes
    self.frame = streamtube(lines, colors=(1,1,1), linewidth=0.1)
    if self.boundary:
      self.scene.add(self.frame)


  def update_atomic_positions ( self ):
    from fury.utils import vertices_from_actor,update_actor

    ainfo,binfo,linfo = self.model.lattice_atoms_bonds(*self.nsc,
                        self.bond_type, self.frame_index,
                        self.constrain_atoms)

    verts = vertices_from_actor(self.atoms)

    diffs = ainfo[0] - self.aposs
    self.aposs = ainfo[0]
    sec = verts.shape[0] // diffs.shape[0]
    for i,d in enumerate(diffs):
      strt,end = i * sec, (i+1) * sec
      verts[strt:end] += d
    update_actor(self.atoms)

    mem = colors_from_actor(self.atoms, 'colors')
    for i,aind in enumerate(self.sel_inds):
      self.set_atom_color(mem, aind, sec, self.sel_cols[i])

    if self.relax_boundary:
      self.scene.rm(self.frame)
      self.update_boundary(*linfo)

    self.update_selections()
    self.scene.ResetCameraClippingRange()
    self.smanager.render()


  def redraw_atomic_model ( self ):
    from fury.utils import update_actor

    self.scene.rm(self.frame)
    if self.atoms is not None:
      self.scene.rm(self.atoms)
    for b in self.bonds:
      self.scene.rm(b)

    self.render_atomic_model()
    self.update_selections()

    self.scene.ResetCameraClippingRange()
    self.smanager.render()


  def update_selections ( self ):

    for b in self.sel_bnds:
      self.scene.rm(b)

    if self.atoms is not None:

      self.sel_forward = True
      sel_inds = self.sel_inds.copy()
      colors = colors_from_actor(self.atoms, 'colors')
      nvert = int(vertices_from_actor(self.atoms).shape[0]/self.natoms)

      self.sel_inds = []
      self.sel_cols = []
      self.sel_bnds = []
      for ind in sel_inds:
        self.selection_logic(colors, ind, nvert)


  def update_atomic_model ( self ):
    if self.model.bonds:
      self.redraw_atomic_model()
    else:
      self.update_atomic_positions()


  def render_atomic_model ( self ):
    from .cylinder import cylinder
    from fury import actor
    import numpy as np

    ainfo,binfo,linfo = self.model.lattice_atoms_bonds(*self.nsc,
                        self.bond_type, self.frame_index,
                        self.constrain_atoms)
 
    self.aposs = ainfo[0]
    self.natoms = ainfo[0].shape[0]

    phi = theta = 30
    if self.natoms > 4000:
      phi = 8; theta = 6
    elif self.natoms > 2000:
      phi = 12; theta = 10
    elif self.natoms > 100:
      phi = theta = 16

    if self.natoms > 0:
      try:
        self.atoms = actor.sphere(ainfo[0], ainfo[1], phi=phi,
                                theta=theta, radii=ainfo[2],
                                use_primitive=False)
      except:
        self.atoms = actor.sphere(ainfo[0], ainfo[1], phi=phi,
                                  theta=theta, radii=ainfo[2])
      self.scene.add(self.atoms)
      self.atoms.AddObserver('LeftButtonPressEvent', self.pick_atom, 1)

    self.bonds = []
    for i in range(binfo[0].shape[0]):
      tbond = cylinder(binfo[0][i], binfo[1][i], binfo[2][i],
              radius=binfo[3][i], heights=binfo[4][i], resolution=20)
      self.scene.add(tbond)
      self.bonds.append(tbond)

    self.update_boundary(*linfo)


  def render_iso_surface ( self, data, origin=(0,0,0), arrows=None, iso_vals=0, colors=(255,110,0), arrow_colors=(255,100,0), arrow_scale=0.25, arrow_anchor='mid', arrow_spacing=0.01, disp_all=False, clip_planes=None, clip_boundary=True):
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
    nsc = self.nsc
    origin = np.array(origin, dtype=float)
    origin -= np.array([((nsc[i]+1)%4)/4 for i in range(3)])
    origin += .5/np.array(data.shape)/nsc*[3-nsc[i]%2 for i in range(3)]
    super().render_iso_surface(self.model.lattice, origin, data, arrows, iso_vals, colors, arrow_colors, arrow_scale, arrow_anchor, arrow_spacing, disp_all, clip_planes, clip_boundary, nsc)


  def start_crystal_view ( self, camera_pos=None, camera_focal=None, camera_up=None ):
    '''
      Begin the render sequence and allow interaction

      Arguments:
        camera_pos (list): 3-vector position to place the camera.
                           None defaults to position central to the
                           rendered objects in the xy plane.
        camera_focal (list): 3-vector positions of the cameras focal
                             point. None defaults to the center of the
                             rendered objects.
        camera_up (list): 3-vector for the cameras "up" orientation.
                          None defaults to unit vector y [0,1,0].
    '''
    super().start_crystal_view(camera_pos, camera_focal, camera_up)
    super().camera_default_position()
    self.smanager.render()
    self.smanager.start()
