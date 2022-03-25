from fury.lib import numpy_support,Actor,CellArray,DataObject,DoubleArray,ImageData,Points,PolyData,PolyDataMapper,Transform,VTK_UNSIGNED_CHAR
from fury.utils import numpy_to_vtk_matrix,set_polydata_colors
import numpy as np

import vtkmodules.vtkCommonCore as ccvtk
import vtkmodules.vtkFiltersCore as fcvtk
import vtkmodules.vtkFiltersGeneral as fgvtk
import vtkmodules.vtkFiltersModeling as fmvtk
import vtkmodules.vtkFiltersGeometry as fgmvtk
import vtkmodules.vtkCommonDataModel as cdmvtk

Planes = cdmvtk.vtkPlanes
Triangle = cdmvtk.vtkTriangle
ClipPolyData = fcvtk.vtkClipPolyData
ContourFilter = fcvtk.vtkContourFilter
MultiThreshold = fgvtk.vtkMultiThreshold
TransformFilter = fgvtk.vtkTransformFilter
UnsignedCharArray = ccvtk.vtkUnsignedCharArray
SelectEnclosedPoints = fmvtk.vtkSelectEnclosedPoints
DataSetSurfaceFilter = fgmvtk.vtkDataSetSurfaceFilter
TransformPolyDataFilter = fgvtk.vtkTransformPolyDataFilter
WindowedSincPolyDataFilter = fcvtk.vtkWindowedSincPolyDataFilter

# Create my iso_surface routine, until I'm allowed to merge the branch into fury
def iso_surface(data, dx, iso_val, origin, colors, bound_planes=None, skew=None):

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
      colors = colors[:-1,:-1,:-1]
      colors = np.ascontiguousarray(np.swapaxes(colors.astype('uint8'),0,2).reshape((np.prod(colors.shape[:3]),3)))
      vtk_colors = numpy_support.numpy_to_vtk(colors, deep=True, array_type=vtk_type)
      vtk_colors.SetNumberOfComponents(3)

    im = ImageData()
    im.SetDimensions(*dims)
    im.SetOrigin(*origin)
    im.SetSpacing(*dx)
    im.AllocateScalars(vtk_type, 1)
    im.GetPointData().SetScalars(uchar_array)

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

    if not one_color:
      im.GetCellData().SetScalars(vtk_colors)

    if bound_planes is None:
        iso_map = PolyDataMapper()
        iso_map.SetInputConnection(iso_transformed.GetOutputPort())

    else:
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

        iso_map = PolyDataMapper()
        iso_map.SetInputConnection(clipper.GetOutputPort())

    if not one_color:
        iso_map.SetScalarModeToUseCellData()
        iso_map.ScalarVisibilityOn()
    else:
        iso_map.ScalarVisibilityOff()

    actor = Actor()
    actor.SetMapper(iso_map)
    if one_color:
        actor.GetProperty().SetColor(*colors/255.)

    return actor
