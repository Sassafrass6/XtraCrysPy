from .XtraCrysPy import XtraCrysPy
from .r_diagram import helix,r_diagram
import numpy as np

class Richardson ( XtraCrysPy ):

  def __init__ ( self, size=(1024, 1024), axes=True, boundary=True, background=(0,0,0), perspective=False, image_prefix='XCP_Image' ):
    super().__init__(size, axes, boundary, background, perspective, image_prefix)

    # Build fake spiral
    nV = 256
    rS = 2
    nCyc = 3
    h = 10
    points = []
    normals = []
    for i in range(nV):
      x = rS * np.cos(2*np.pi*nCyc*i/(nV-1))
      y = rS * np.sin(2*np.pi*nCyc*i/(nV-1))
      z = h * i / nV
      p = np.array([x,y,z])
      points.append(p)
      normals.append(np.array([0,0,z])-p)

    # UI SETUP

    # Build spirals like this...
    #r_actor = r_diagram(points, normals, .5, [1,.3,0])
    #self.scene.add(r_actor)

    # Build a helix like this?
    h_actor = helix([[-1,-2,-3],[3,2,1]], 4, 3, .5, [1,.3,0])
    self.scene.add(h_actor)

    # Finish setting up the scene
    self.scene.ResetCamera()
    self.smanager.render()

    # User should call this
    self.start_crystal_view()
