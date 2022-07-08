from fury.lib import (numpy_support, PolyData, Points, Actor,
                      CylinderSource, PolyDataMapper, Glyph3D,
                      VTK_DOUBLE, Transform)
from fury.utils import (set_polydata_vertices, set_polydata_triangles,
                        numpy_to_vtk_matrix, numpy_to_vtk_colors,
                        numpy_to_vtk_points)
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter

import numpy as np

def repeat_sources(centers, colors, active_scalars=1., directions=None,
                   source=None, vertices=None, faces=None, orientation=None):
  """Transform a vtksource to glyph."""
  if source is None and faces is None:
    raise IOError("A source or faces should be defined")

  if np.array(colors).ndim == 1:
    colors = np.tile(colors, (len(centers), 1))

  pts = numpy_to_vtk_points(np.ascontiguousarray(centers))
  cols = numpy_to_vtk_colors(255 * np.ascontiguousarray(colors))
  cols.SetName('colors')
  if isinstance(active_scalars, (float, int)):
    active_scalars = np.tile(active_scalars, (len(centers), 1))
  if isinstance(active_scalars, np.ndarray):
    ascalars = numpy_support.numpy_to_vtk(np.asarray(active_scalars),
                                          deep=True,
                                          array_type=VTK_DOUBLE)
    ascalars.SetName('active_scalars')

  if directions is not None:
    directions_fa = numpy_support.numpy_to_vtk(np.asarray(directions),
                                               deep=True,
                                               array_type=VTK_DOUBLE)
    directions_fa.SetName('directions')

  polydata_centers = PolyData()
  polydata_geom = PolyData()

  if faces is not None:
    directions_fa = numpy_support.numpy_to_vtk(np.asarray(directions),
                                               deep=True,
                                               array_type=VTK_DOUBLE)
    directions_fa.SetName('directions')

  polydata_centers = PolyData()
  polydata_geom = PolyData()

  if faces is not None:
    set_polydata_vertices(polydata_geom, vertices)
    set_polydata_triangles(polydata_geom, faces)

  polydata_centers.SetPoints(pts)
  polydata_centers.GetPointData().AddArray(cols)
  if directions is not None:
    polydata_centers.GetPointData().AddArray(directions_fa)
    polydata_centers.GetPointData().SetActiveVectors('directions')
  if isinstance(active_scalars, np.ndarray):
    polydata_centers.GetPointData().AddArray(ascalars)
    polydata_centers.GetPointData().SetActiveScalars('active_scalars')

  glyph = Glyph3D()
  if faces is None:
    if orientation is not None:
      transform = Transform()
      transform.SetMatrix(numpy_to_vtk_matrix(orientation))
      rtrans = vtkTransformPolyDataFilter()
      rtrans.SetInputConnection(source.GetOutputPort())
      rtrans.SetTransform(transform)
      source = rtrans
      glyph.SetSourceConnection(source.GetOutputPort())
    else:
      glyph.SetSourceData(polydata_geom)
  glyph.SetInputData(polydata_centers)
  glyph.SetOrient(True)
  glyph.SetScaleModeToScaleByScalar()
  glyph.SetVectorModeToUseVector()
  glyph.Update()

  mapper = PolyDataMapper()
  mapper.SetInputData(glyph.GetOutput())
  mapper.SetScalarModeToUsePointFieldData()
  mapper.SelectColorArray('colors')

  actor = Actor()
  actor.SetMapper(mapper)
  return actor

def cylinder(centers, directions, colors, radius=0.05, heights=1,
             capped=False, resolution=6, vertices=None, faces=None):

  if faces is None:
    src = CylinderSource()
    src.SetCapping(capped)
    src.SetResolution(resolution)
    src.SetRadius(radius)
    rotate = np.array([[0, 1, 0, 0],
                       [-1, 0, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]])
  else:
    src = None
    rotate = None

  cylinder_actor = repeat_sources(centers=centers, colors=colors,
                                    directions=directions,
                                    active_scalars=heights, source=src,
                                    vertices=vertices, faces=faces,
                                    orientation=rotate)

  return cylinder_actor

