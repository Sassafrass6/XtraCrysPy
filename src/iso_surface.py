from fury.lib import numpy_support,Actor,CellArray,DataObject,ImageData,Points,PolyData,PolyDataMapper,Transform,VTK_UNSIGNED_CHAR
from fury.utils import numpy_to_vtk_matrix,set_polydata_colors
import numpy as np

import vtkmodules.vtkCommonCore as ccvtk
import vtkmodules.vtkFiltersCore as fcvtk
import vtkmodules.vtkFiltersGeneral as fgvtk
import vtkmodules.vtkFiltersModeling as fmvtk
import vtkmodules.vtkFiltersGeometry as fgmvtk
import vtkmodules.vtkCommonDataModel as cdmvtk

Triangle = cdmvtk.vtkTriangle
ContourFilter = fcvtk.vtkContourFilter
FlyingEdges3D = fcvtk.vtkFlyingEdges3D
MultiThreshold = fgvtk.vtkMultiThreshold
TransformFilter = fgvtk.vtkTransformFilter
UnsignedCharArray = ccvtk.vtkUnsignedCharArray
SelectEnclosedPoints = fmvtk.vtkSelectEnclosedPoints
DataSetSurfaceFilter = fgmvtk.vtkDataSetSurfaceFilter
TransformPolyDataFilter = fgvtk.vtkTransformPolyDataFilter
WindowedSincPolyDataFilter = fcvtk.vtkWindowedSincPolyDataFilter

# Create my iso_surface routine, until I'm allowed to merge the branch into fury
def iso_surface(data, iso_val, origin, colors, bound_polys=None, lattice=None):

    if data.ndim != 3:
        raise ValueError('Only 3D arrays are currently supported.')

    dims = data.shape
    npnt = np.prod(dims)
    voxsz = [2/s for s in dims]

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
    im.SetSpacing(voxsz[0], voxsz[1], voxsz[2])
    im.AllocateScalars(vtk_type, 1)
    im.GetPointData().SetScalars(uchar_array)

    iso = ContourFilter()
    iso.SetInputData(im)
    iso.SetValue(0, iso_val)

    lpad = np.eye(4)
    lpad[:3, :3] = lattice.T
    transform = Transform()
    transform.SetMatrix(numpy_to_vtk_matrix(lpad))

    latT = TransformPolyDataFilter()
    latT.SetInputConnection(iso.GetOutputPort())
    latT.SetTransform(transform)

    if not one_color:
      im.GetCellData().SetScalars(vtk_colors)

    if bound_polys is None:
        iso_map = PolyDataMapper()
        iso_map.SetInputConnection(latT.GetOutputPort())

    else:
        points = Points()
        for p in bound_polys[0]:
            points.InsertNextPoint(*p)

        tris = CellArray()
        for vtx in bound_polys[1]:
            t = Triangle()
            for i,v in enumerate(vtx):
                t.GetPointIds().SetId(i, v)
            tris.InsertNextCell(t)

        poly = PolyData()
        poly.SetPoints(points)
        poly.SetPolys(tris)

        select = SelectEnclosedPoints()
        select.SetInputConnection(latT.GetOutputPort())
        select.SetSurfaceData(poly)

        thresh = MultiThreshold()
        inside = thresh.AddBandpassIntervalSet(1, 1, DataObject.FIELD_ASSOCIATION_POINTS, 'SelectedPoints', 0, 1)
        thresh.SetInputConnection(select.GetOutputPort())
        thresh.OutputSet(inside)
        thresh.Update()

        surface = DataSetSurfaceFilter()
        surface.SetInputData(thresh.GetOutput().GetBlock(inside).GetBlock(0))

        smooth = WindowedSincPolyDataFilter()
        smooth.SetInputConnection(surface.GetOutputPort())
        smooth.SetNumberOfIterations(30)
        smooth.BoundarySmoothingOn()
        smooth.SetEdgeAngle(180)
        smooth.GetOutput().GetCellData().SetScalars(iso.GetOutput().GetCellData().GetScalars())
        smooth.Update()

        iso_map = PolyDataMapper()
        iso_map.SetInputConnection(smooth.GetOutputPort())

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
