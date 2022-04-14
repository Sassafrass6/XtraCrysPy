
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

