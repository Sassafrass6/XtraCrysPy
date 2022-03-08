from fury.lib import numpy_support,Actor,CellArray,DataObject,DataSetMapper,ImageData,Points,PolyData,PolyDataMapper
import numpy as np

import vtkmodules.vtkFiltersCore as fcvtk
import vtkmodules.vtkFiltersGeneral as fgvtk
import vtkmodules.vtkFiltersModeling as fmvtk
import vtkmodules.vtkCommonDataModel as cdmvtk

Triangle = cdmvtk.vtkTriangle
FlyingEdges3D = fcvtk.vtkFlyingEdges3D
MultiThreshold = fgvtk.vtkMultiThreshold
SelectEnclosedPoints = fmvtk.vtkSelectEnclosedPoints

# Create my iso_surface routine, until I'm allowed to merge the branch into fury
def iso_surface(data, iso_val, origin, color, bound_polys=None):

    if data.ndim != 3:
        raise ValueError('Only 3D arrays are currently supported.')

    nb_components = 1
    dims = data.shape
    voxsz = [1/s for s in dims]

    data = data.astype('uint8')
    vtk_type = numpy_support.get_vtk_array_type(data.dtype)
    data = np.ascontiguousarray(np.swapaxes(data, 0, 2)).ravel()
    uchar_array = numpy_support.numpy_to_vtk(data, deep=0)

    im = ImageData()
    im.SetDimensions(*dims)
    im.SetOrigin(*origin)
    im.SetSpacing(voxsz[0], voxsz[1], voxsz[2])
    im.AllocateScalars(vtk_type, nb_components)
    im.GetPointData().SetScalars(uchar_array)

    iso = FlyingEdges3D()
    iso.SetInputData(im)
    iso.SetValue(0, iso_val)

    if bound_polys is None:
        iso_map = PolyDataMapper()
        iso_map.SetInputConnection(iso.GetOutputPort())
        iso_map.ScalarVisibilityOff()

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
        select.SetInputConnection(iso.GetOutputPort())
        select.SetSurfaceData(poly)

        thresh = MultiThreshold()
        inside = thresh.AddBandpassIntervalSet(1, 1, DataObject.FIELD_ASSOCIATION_POINTS, 'SelectedPoints', 0, 1)
        thresh.SetInputConnection(select.GetOutputPort())
        thresh.OutputSet(inside)
        thresh.Update()

        iso_map = DataSetMapper()
        iso_map.SetInputData(thresh.GetOutput().GetBlock(inside).GetBlock(0))
        iso_map.ScalarVisibilityOff()

    actor = Actor()
    actor.SetMapper(iso_map)
    actor.GetProperty().SetColor(*color)

    return actor
