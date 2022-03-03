import numpy as np

class XtraCrysPy:

  def __init__ ( self, size=(1024, 1024), structure=None ):
    '''

    Arguments:
    '''
    from fury import window

    self.scene = window.Scene()
    self.smanager = window.ShowManager(self.scene, size=size, order_transparent=True)
    self.smanager.initialize()

    self.axes = None

  def start_crystal_view ( self ):

    self.smanager.start()


  def render_atomic_model ( self, model ):
    '''
    '''
    from fury import actor
    import numpy as np

    ainfo,binfo = model.lattice_atoms_bonds(1,1,1)

    if ainfo[0].shape[0] > 0:
      self.atoms = actor.sphere(centers=ainfo[0], colors=ainfo[1], radii=ainfo[2])
      self.scene.add(self.atoms)

    self.bonds = []
    for i in range(0, binfo[0].shape[0], 2):
      tbond = actor.cylinder(binfo[0][i], binfo[1][i], binfo[2][i], radius=binfo[3][i], heights=binfo[4][i], resolution=10)
      self.scene.add(tbond)
      self.bonds.append(tbond)
    '''
    xyz = np.array([[2,0,0], [0,2,0]])
    colors = np.array([[1,0,0,1], [1,1,0,1]])
    radii = np.array([.5,.5])
    self.atoms = actor.sphere(centers=xyz, colors=colors, radii=radii)
    self.scene.add(self.atoms)

    con = xyz[1]-xyz[0]
    cen = (xyz[0]+xyz[1])/2
    bpos = np.array([cen-con/4, cen+con/4])
    height = np.linalg.norm(con)/2
    bdir = np.cross(con,np.cross(xyz[0],xyz[1]))
    bdir = np.array([bdir, bdir])
    color = colors
    radius = (radii[0]+radii[1])/8
    self.bonds = actor.cylinder(bpos, bdir, color, radius=radius, heights=height)
    self.scene.add(self.bonds)
    '''


  def show_axes ( self, show=True ):
    '''
    '''
    return
    from fury import actor

    self.axes = actor.axes()
    self.scene.add(self.axes)
    self.scene.reset_camera()
