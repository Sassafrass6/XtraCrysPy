from fury.utils import colors_from_actor, update_actor
from fury.utils import vertices_from_actor
from .XtraCrysPy import XtraCrysPy
import numpy as np

class XCP_Atoms ( XtraCrysPy ):

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True,
                 background=(0,0,0), perspective=False, model=None,
                 params={}, relax=False, nsc=(1,1,1),
		 bond_type='Stick', sel_type='Chain', unit='angstrom',
		 runit='degree' ):
    super().__init__(size, axes, boundary, background, perspective)
    from .Model import Model

    self.nsc = nsc
    self.sel_forward = True
    self.sel_type = sel_type
    self.bond_type = bond_type
    self.constrain_atoms = True

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
        model = Model(params, fname=model, relax=relax)
      else:
        s = 'Argument \'model\' must be a file name or a Model object.'
        raise TypeError(s)
    elif params:
      model = Model(params, fname=None, relax=relax)
    else:
      raise Exception('Must specify model or params arguments')

    self.model = model
    self.relax = relax
    self.relax_index = 0
    self.nrelax = self.model.atoms.shape[0]

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


  def setup_ui ( self ):
    from fury.data import read_viz_icons
    from fury import ui

    size = self.wsize
    icons = [('left',read_viz_icons(fname='circle-left.png')), 
              ('down',read_viz_icons(fname='circle-down.png')),
              ('right',read_viz_icons(fname='circle-right.png'))]

    self.sel_type_button = ui.Button2D(icon_fnames=icons[2:0:-1],
                                       size=(40,40))

    sel_types = ['Chain', 'Info', 'Angle', 'Distance']
    if self.sel_type not in sel_types:
      s = '''{} is not a valid selection type. Choose from:
               Info, Angle Chain, or Distance.'''
      self.sel_type = 'Chain'
      print(s)

    self.sel_type_menu = ui.ListBox2D(sel_types, multiselection=False,
                         font_size=24, line_spacing=2, size=(150,220),
                         background_opacity=.9)

    sel = sel_types.index(self.sel_type)
    self.sel_type_menu.select(self.sel_type_menu.slots[sel])
    self.sel_type_menu.on_change = self.update_selection_type

    self.sel_menu_vis = False
    sel_button = self.sel_type_button
    self.sel_type_menu.set_visibility(False)
    sel_button.on_left_mouse_button_clicked = self.toggle_sel_menu

    checkbox = ['Constrain']
    self.constrain_checkbox = ui.Checkbox(checkbox, checkbox,
                              font_size=24, font_family='Arial',
                              position=(10,size[1]-65))
    self.scene.add(self.constrain_checkbox)
    self.constrain_checkbox.on_change = self.update_constrain

    self.sel_panel = ui.Panel2D((40,40), (10, size[1]-110), opacity=0)
    self.sel_panel.add_element(self.sel_type_button, (0,0))
    self.sel_panel.add_element(self.sel_type_menu, (50,-210))
    self.scene.add(self.sel_panel)

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

      self.relax_panel = ui.Panel2D((110,100),
                         (size[0]-120,size[1]-60), (0,0,0))
      self.relax_panel.add_element(left_button, (0,0))
      self.relax_panel.add_element(right_button, (60,0))

      stext = '1/{}'.format(self.nrelax)
      self.relax_text = ui.TextBlock2D(text=stext, 
                        justification='center', font_size=30,
                        vertical_justification='top')
      self.relax_panel.add_element(self.relax_text, (55,-10))
      self.scene.add(self.relax_panel)


  def update_buttons ( self, caller, event ):
    super().update_buttons(caller, event)
    x,y = self.scene.GetSize()
    self.sel_tpanel.center = (x-160, 50)
    self.sel_panel.position = (10, y-110)
    self.constrain_checkbox.position = (10, y-65)
    if self.relax:
      self.relax_panel.position = (x-120, y-60)


  def update_selection_type ( self ):
    self.clear_selection_text()
    self.sel_type = self.sel_type_menu.selected[0] 
    for _ in range(len(self.sel_bnds)):
      self.pop_sbond()
    colors = colors_from_actor(self.atoms)
    nvert = int(vertices_from_actor(self.atoms).shape[0]/self.natoms)
    for i in self.sel_inds.copy():
      self.pop_atom(colors, i, nvert)
    update_actor(self.atoms)


  def clear_selection_text ( self ):
    self.sel_text.message = ''
    self.sel_tpanel.set_visibility(False)


  def update_relax_text ( self ):
    text = '{}/{}'.format(self.relax_index+1, self.nrelax)
    self.relax_text.message = text


  def update_selection_text ( self, text ):
    self.sel_text.message = text
    self.sel_tpanel.set_visibility(True)


  def update_constrain ( self, checkboxes ):
    self.constrain_atoms = not self.constrain_atoms
    self.update_atomic_model()


  def toggle_sel_menu ( self, iren, caller, event ):
    self.sel_menu_vis = not self.sel_menu_vis
    self.sel_type_menu.set_visibility(self.sel_menu_vis)
    self.sel_type_button.next_icon()
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
      brad = 0.01 + brad / 2
      tbond = actor.cylinder([cent], [conn], [self.scolor/255],
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
        spec = self.model.species[index%self.model.natoms]
        message = 'Atom {} : {}'.format(index, spec)
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


  def relax_forward ( self, iren, obj, event ):
    if self.relax_index < self.nrelax-1:
      self.relax_index += 1
      self.update_relax_text()
      self.update_atomic_model()


  def relax_backward ( self, iren, obj, event ):
    if self.relax_index > 0:
      self.relax_index -= 1
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


  def update_atomic_model ( self ):
    from fury.utils import update_actor

    self.scene.rm(self.frame)
    self.scene.rm(self.atoms)
    for b in self.bonds:
      self.scene.rm(b)
    for b in self.sel_bnds:
      self.scene.rm(b)

    self.render_atomic_model()

    self.sel_forward = True
    sel_inds = self.sel_inds.copy()
    colors = colors_from_actor(self.atoms, 'colors')
    nvert = int(vertices_from_actor(self.atoms).shape[0]/self.natoms)

    self.sel_inds = []
    self.sel_cols = []
    self.sel_bnds = []
    for ind in sel_inds:
      self.selection_logic(colors, ind, nvert)

    self.smanager.render()


  def render_atomic_model ( self ):
    '''
    '''
    from fury import actor
    import numpy as np

    ainfo,binfo,linfo = self.model.lattice_atoms_bonds(*self.nsc,
                        self.bond_type, self.relax_index,
                        self.constrain_atoms)
 
    self.aposs = ainfo[0]
    self.natoms = ainfo[0].shape[0]
    if self.natoms > 0:
      self.atoms = actor.sphere(centers=ainfo[0],
                                colors=ainfo[1], radii=ainfo[2])
      self.scene.add(self.atoms)
      self.atoms.AddObserver('LeftButtonPressEvent', self.pick_atom, 1)

    self.bonds = []
    for i in range(binfo[0].shape[0]):
      tbond = actor.cylinder(binfo[0][i], binfo[1][i], binfo[2][i],
              radius=binfo[3][i], heights=binfo[4][i], resolution=20)
      self.scene.add(tbond)
      self.bonds.append(tbond)

    self.bound_planes = linfo[0]
    self.frame = actor.streamtube(linfo[1], colors=(1,1,1), linewidth=0.1)
    if 'Boundary' in self.frame_checkbox.checked_labels:
      self.scene.add(self.frame)

    self.scene.ResetCamera()


  def render_iso_surface ( self, data, iso_vals=0, colors=(255,110,0), disp_all=False ):
    nsc = self.nsc
    origin = -np.array([((nsc[i]+1)%4)/4 for i in range(3)])
    origin += .5/np.array(data.shape)/nsc*[3-nsc[i]%2 for i in range(3)]
    super().render_iso_surface(self.model.lattice, origin, data, iso_vals, colors, disp_all, nsc)
