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
ContourFilter = fcvtk.vtkContourFilter
MultiThreshold = fgvtk.vtkMultiThreshold
TransformFilter = fgvtk.vtkTransformFilter
UnsignedCharArray = ccvtk.vtkUnsignedCharArray
SelectEnclosedPoints = fmvtk.vtkSelectEnclosedPoints
DataSetSurfaceFilter = fgmvtk.vtkDataSetSurfaceFilter
TransformPolyDataFilter = fgvtk.vtkTransformPolyDataFilter
WindowedSincPolyDataFilter = fcvtk.vtkWindowedSincPolyDataFilter

def iso_surface(data, dx, iso_val, origin, colors, bound_planes=None, skew=None, arrows=None):

    if data.ndim != 3:
        raise ValueError('Only 3D arrays are currently supported.')

    dims = data.shape
    npnt = np.prod(dims)

    data = data.astype('uint8')
    vtk_type = numpy_support.get_vtk_array_type(data.dtype)

    data = np.ascontiguousarray(np.swapaxes(data,0,2).reshape(npnt))
    uchar_array = numpy_support.numpy_to_vtk(data, deep=False, array_type=vtk_type)

    one_color = True
    if len(colors.shape) != 1:
      one_color = False

    im = ImageData()
    im.SetDimensions(*dims)
    im.SetOrigin(*origin)
    im.SetSpacing(*dx)
    im.AllocateScalars(vtk_type, 1)
    im.GetPointData().SetScalars(uchar_array)

    if arrows is not None:
      arrows = np.swapaxes(0.01*arrows,0,2).reshape((npnt,3))
      adirs = FloatArray()
      adirs.SetName('directions')
      adirs.SetNumberOfComponents(3)
      adirs.SetNumberOfTuples(npnt)
      for i in range(npnt):
        adirs.SetTuple(i, arrows[i])

      im.GetPointData().AddArray(adirs)
      im.GetPointData().SetActiveVectors('directions')

    if not one_color:
      if arrows is not None:
        tcol = np.ascontiguousarray(np.swapaxes(colors.astype('uint8'),0,2).reshape((np.prod(colors.shape[:3]),3)))
        acolors = numpy_support.numpy_to_vtk(tcol, deep=True, array_type=vtk_type)
        acolors.SetNumberOfComponents(3)
        acolors.SetName('colors')
        im.GetPointData().AddArray(acolors)

      colors = colors[:-1,:-1,:-1]
      colors = np.ascontiguousarray(np.swapaxes(colors.astype('uint8'),0,2).reshape((np.prod(colors.shape[:3]),3)))
      vtk_colors = numpy_support.numpy_to_vtk(colors, deep=True, array_type=vtk_type)
      vtk_colors.SetNumberOfComponents(3)
      im.GetCellData().SetScalars(vtk_colors)

    iso = ContourFilter()
    iso.SetInputData(im)
    iso.SetValue(0, iso_val)

    lpad = np.eye(4)
    lpad[:3, :3] = skew.T
    transform = Transform()
    transform.SetMatrix(numpy_to_vtk_matrix(lpad))

    iso_transformed = TransformPolyDataFilter()
    iso_transformed.SetInputConnection(iso.GetOutputPort())
    iso_transformed.SetTransform(transform)

    if bound_planes is not None:
        clip_centers = Points()
        clip_normals = DoubleArray()
        clip_normals.SetNumberOfComponents(3)
        for i in range(bound_planes.shape[1]):
          clip_centers.InsertNextPoint(*bound_planes[0,i])
          clip_normals.InsertNextTuple3(*bound_planes[1,i])
        clip_planes = Planes()
        clip_planes.SetPoints(clip_centers)
        clip_planes.SetNormals(clip_normals)

        clipper = ClipPolyData()
        clipper.SetInputConnection(iso_transformed.GetOutputPort())
        clipper.SetClipFunction(clip_planes)
        clipper.InsideOutOn()

        iso_transformed = clipper

    iso_map = PolyDataMapper()
    iso_map.SetInputConnection(iso_transformed.GetOutputPort())

    aactor = None
    if arrows is not None:

      arrow = ArrowSource()
      arrow_glyph = Glyph3D()
      arrow_glyph.SetInputData(iso_transformed.GetOutput())
      arrow_glyph.OrientOn()
      arrow_glyph.SetSourceConnection(arrow.GetOutputPort())
      arrow_glyph.SetScaleModeToDataScalingOff()
      arrow_glyph.SetScaleFactor(0.025)
      arrow_glyph.Update()

      amapper = PolyDataMapper()
      amapper.SetInputConnection(arrow_glyph.GetOutputPort())

      aactor = Actor()
      aactor.SetMapper(amapper)

    if not one_color:
        iso_map.SetScalarModeToUseCellData()
        iso_map.ScalarVisibilityOn()

        if arrows is not None:
          amapper.SetScalarModeToUsePointFieldData()
          amapper.SelectColorArray('colors')
          amapper.ScalarVisibilityOn()
    else:
        iso_map.ScalarVisibilityOff()

        if arrows is not None:
          amapper.ScalarVisibilityOff()

    actor = Actor()
    actor.SetMapper(iso_map)
    if one_color:
        actor.GetProperty().SetColor(*colors/255.)

    return actor, aactor

