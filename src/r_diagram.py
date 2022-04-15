
from fury.lib import numpy_support,Actor,CellArray,DataObject,DoubleArray,ImageData,Points,PolyData,PolyDataMapper,Transform,VTK_UNSIGNED_CHAR
from fury.utils import numpy_to_vtk_matrix,set_polydata_colors
import numpy as np

import vtkmodules.vtkCommonCore as ccvtk
import vtkmodules.vtkFiltersCore as fcvtk
import vtkmodules.vtkFiltersGeneral as fgvtk
import vtkmodules.vtkFiltersSources as fsvtk
import vtkmodules.vtkFiltersModeling as fmvtk
import vtkmodules.vtkFiltersGeometry as fgmvtk
import vtkmodules.vtkCommonDataModel as cdmvtk

Planes = cdmvtk.vtkPlanes
Glyph3D = fcvtk.vtkGlyph3D
Triangle = cdmvtk.vtkTriangle
FloatArray = ccvtk.vtkFloatArray
ArrowSource = fsvtk.vtkArrowSource
LookupTable = ccvtk.vtkLookupTable
ClipPolyData = fcvtk.vtkClipPolyData
RibbonFilter = fmvtk.vtkRibbonFilter
ContourFilter = fcvtk.vtkContourFilter
MultiThreshold = fgvtk.vtkMultiThreshold
TransformFilter = fgvtk.vtkTransformFilter
UnsignedCharArray = ccvtk.vtkUnsignedCharArray
SelectEnclosedPoints = fmvtk.vtkSelectEnclosedPoints
DataSetSurfaceFilter = fgmvtk.vtkDataSetSurfaceFilter
TransformPolyDataFilter = fgvtk.vtkTransformPolyDataFilter
WindowedSincPolyDataFilter = fcvtk.vtkWindowedSincPolyDataFilter

def helix ( center_line, ncyc, radius, width, color, npoints=100, pref=(1,0,0)):
  '''
    Arguments:
      center_line (list): (N,3) array of xyz points on the helix central line, where N is the number of points in each ribbon.
      ncyc (float): The number of rotations in the helix
      radius (float): Radius of the helix
      width (float): Ribbon width.
      color (list): RGB color for the ribbon
      npoints (int): Number of points in the helix descretization
  '''

  if ncyc < .5:
    raise ValueError('Must allow more than 1/4 rotation for each endcap (ncyc>.5)')
  center_line = np.array(center_line)
  if center_line.shape[0] != 2:
    raise ValueError('ONLY DOING 2pt LINE NOW')

  p0 = center_line[0]
  p1 = center_line[1]
  pdiff = p1 - p0
  pref = np.array(pref)
  plen = np.linalg.norm(pdiff)
  nx = lambda pt : pdiff / plen # Normal
  px = lambda pt : p0 + pt * pdiff # Point
  tx = lambda pt : pt * 2 * np.pi * ncyc # Phase
  if np.isclose(plen, 0):
    raise ValueError('Line contains identical points')

  def perpendicular ( v1, v2 ):
    perp = np.cross(v1, v2)
    return perp / np.linalg.norm(perp)

  width /= 2
  hwidths = []
  hpoints = []
  hnormals = []

  hstart = 1 / ncyc / 4
  hend = 1 - hstart
  npcaps = int(npoints * hstart)
  for p in np.linspace(0, hstart, npcaps*10)[:-1]:
    phase = tx(p)
    center = px(p)
    normal = nx(p)
    perp = perpendicular(pref, normal)
    vrot = perp*np.cos(phase) + np.cross(normal, perp) * np.sin(phase)
    vrot += normal * normal.dot(perp) * (1 - np.cos(phase))
    nn = -vrot
    hp = center + vrot * np.sin(phase)
    hpoints.append(hp)
    hnormals.append(nn)
    hwidths.append(p*width/hstart)

  for p in np.linspace(hstart, hend, npoints - 2*npcaps):
    phase = tx(p)
    center = px(p)
    normal = nx(p)
    perp = perpendicular(pref, normal)
    vrot = perp*np.cos(phase) + np.cross(normal, perp) * np.sin(phase)
    vrot += normal * normal.dot(perp) * (1 - np.cos(phase))

    nn = -vrot
    hp = center + vrot
    hpoints.append(hp)
    hnormals.append(nn)
    hwidths.append(width)

  width = hwidths
  npoints = len(hpoints)
  r_pos = Points()
  r_lines = CellArray()
  r_lines.InsertNextCell(npoints)
  for i,p in enumerate(hpoints):
    r_pos.InsertPoint(i, *p)
    r_lines.InsertCellPoint(i)

  r_poly = PolyData()
  r_poly.SetPoints(r_pos)
  r_poly.SetLines(r_lines)

  r_norms = FloatArray()
  r_norms.SetNumberOfComponents(3)
  r_norms.SetNumberOfTuples(len(hnormals))
  for i,n in enumerate(hnormals):
    r_norms.SetTuple(i, n)
  r_poly.GetPointData().SetNormals(r_norms)

  vary_width = not isinstance(width, (int,float))
  if vary_width:
    r_wid = FloatArray()
    r_wid.SetNumberOfComponents(1)
    r_wid.SetNumberOfTuples(len(hpoints))
    for i,w in enumerate(width):
      r_wid.SetTuple(i, [w])
    r_poly.GetPointData().SetScalars(r_wid)

  r_ribbon = RibbonFilter()
  r_ribbon.SetInputData(r_poly)
  r_ribbon.UseDefaultNormalOff()
  if vary_width:
    r_ribbon.VaryWidthOn()
  else:
    r_ribbon.VaryWidthOff()
    r_ribbon.SetWidth(width)

  r_map = PolyDataMapper()
  r_map.SetInputConnection(r_ribbon.GetOutputPort())
  r_map.ScalarVisibilityOff()

  r_actor = Actor()
  r_actor.SetMapper(r_map)
  r_actor.GetProperty().SetColor(*color)

  return r_actor


def r_diagram(points, normals, width, color):
  '''
    Arguments:
      points (list): (N,3) array of points for the ribbon, N is the number of points in each ribbon.
      width (float): Ribbon width.
      color (list): RGB color for the ribbon
      !!!!colors (list): (N,3) or (N,4) list to specify RGB or RGBA colors
  '''

  nlines = len(points)

  r_pos = Points()
  r_lines = CellArray()
  r_lines.InsertNextCell(nlines)
  for i,p in enumerate(points):
    r_pos.InsertPoint(i, *p)
    r_lines.InsertCellPoint(i)

  r_poly = PolyData()
  r_poly.SetPoints(r_pos)
  r_poly.SetLines(r_lines)

  r_norms = FloatArray()
  r_norms.SetNumberOfComponents(3)
  r_norms.SetNumberOfTuples(len(normals))
  for i,n in enumerate(normals):
    r_norms.SetTuple(i, n)
  r_poly.GetPointData().SetNormals(r_norms)

  vary_width = not isinstance(width, (int,float))
  if vary_width:
    r_wid = FloatArray()
    r_wid.SetNumberOfComponents(1)
    r_wid.SetNumberOfTuples(len(points))
    for i,w in enumerate(width):
      r_wid.SetTuple(i, [w])
    r_poly.GetPointData().SetScalars(r_wid)

  r_ribbon = RibbonFilter()
  r_ribbon.SetInputData(r_poly)
  r_ribbon.UseDefaultNormalOff()
  if vary_width:
    r_ribbon.VaryWidthOn()
  else:
    r_ribbon.VaryWidthOff()
    r_ribbon.SetWidth(width)

  r_map = PolyDataMapper()
  r_map.SetInputConnection(r_ribbon.GetOutputPort())
  r_map.ScalarVisibilityOff()

  r_actor = Actor()
  r_actor.SetMapper(r_map)
  r_actor.GetProperty().SetColor(*color)

  return r_actor

