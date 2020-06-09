bl_info = {
    "name": "xcp reader",
    "blender": (2, 83, 0),
    "category": "Import-Export",
}

"""A WIP plugin for blender for reading xcp objects (XML-like)
Current design uses dictionaries over objects for now, since this is
a more direct interpretation of XML
"""

import bpy
import numpy as np # I think this comes with blender by default, doublecheck this
from xcp_io import read_xcp

def add_camera(position, rotation):
    cam = bpy.data.cameras.new("Camera 1")
    cam.lens = 18

    # create the first camera object
    objref = bpy.data.objects.new("Camera 1", cam1)
    objref.location = position
    objref.rotation_euler = rotation
    return objref

def draw_ball(position, radius):
    # TODO add this to a collection of sorts
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=64, ring_count=48, 
        align='WORLD', enter_editmode=False, 
        location=(0., 0., 0.), rotation=(0., 0., 0.))
        
    ball = bpy.context.view_layer.objects.active
    ball.location = position
    ball.scale = radius
    return ball

def draw_stick(origin, vec1):
    # next, get rotation matrix
    norm = np.linalg.norm(vec1)
    vec1 /= norm
    vec2 = np.array([0, 0, 1])
    vec3 = np.cross(vec1, vec2)
    dot = np.dot(vec1, vec2)
    skew = np.array([[0, -vec3[2], vec3[1]],
                     [vec3[2], 0, -vec3[0]],
                     [-vec3[1], vec3[0], 0]])
    rotator = np.identity(3) + skew + (np.matmul(skew, skew) * (1.0 / (1.0 + dot)))

    # convert to 4-matrix
    rot4 = np.zeros((4, 4), dtype=float)
    rot4[0:3, 0:3] = rotator
    rot4[3, 3] = 1

    # now use this rotation matrix to rotate the bond
    bpy.ops.mesh.primitive_cylinder_add(location=tuple(origin))
    stick = bpy.context.view_layer.objects.active

    # scale accordingly
    scale4 = np.eye(4, dtype=float)
    scale4[0, 0] = 0.2
    scale4[1, 1] = 0.2
    scale4[2, 2] = norm / 2

    # apply these transformations
    stick.data.transform(scale4)
    stick.data.transform(rot4)
    stick.data.update()

    return stick

def draw_atom(atom):
    objref = draw_ball(atom["pos"], atom["spinfo"]["scale"])
    return objref

def draw_bond(bond):
    pointA = np.array(bond["A"]["pos"])
    pointB = np.array(bond["B"]["pos"])
    objref = []

    midpoint = (pointA + pointB) / 2.0
    originA = (pointA + midpoint) / 2.0
    originB = (midpoint + pointB) / 2.0
    vecA = pointA - midpoint
    vecB = midpoint - pointB
    objref.append(draw_stick(originA, vecA))
    objref.append(draw_stick(originB, vecB))
    return objref

def create_material(id, color):
    mat = bpy.data.materials.new('Material.{}'.format(id))
    mat.diffuse_color = color
    mat.specular_intensity = 0
    return mat

def init_materials(species):
    for key in species:
        spec = species[key]
        spec["material"] = create_material(key, spec["color"])

def import_xcp(context, filepath, use_some_setting):
    xcp = read_xcp(filepath)

    # do some processing in between perhaps and then finally here, draw every atom
    init_materials(xcp["SPECIES"])
    for key in xcp["ATOMS"]:
        atom = xcp["ATOMS"][key]
        atomobj = draw_atom(atom)
        # link materials (expects init_materials to have been called)
        atomobj.data.materials.append(atom["spinfo"]["material"])
    for key in xcp["BONDS"]:
        bond = xcp["BONDS"][key]
        bondobj = draw_bond(bond)
        bondobj[0].data.materials.append(bond["A"]["spinfo"]["material"])
        bondobj[1].data.materials.append(bond["B"]["spinfo"]["material"])

    if "CAMERA" in xcp["SCENE"]:
        camera = xcp["SCENE"]["CAMERA"]
        cameraobj = add_camera(camera["pos"], camera["rot"])
        
    # tells blender that the addon main function is successful
    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator


class ImportXCP(Operator, ImportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "importx.xcp"
    bl_label = "Import XCP Data"

    # ImportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob: StringProperty(
        default="*.xml",
        options={'HIDDEN'},
        maxlen=255,  # Max internal buffer length, longer would be clamped.
    )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    use_setting: BoolProperty(
        name="Example Boolean",
        description="Example Tooltip",
        default=True,
    )

    type: EnumProperty(
        name="Example Enum",
        description="Choose between two items",
        items=(
            ('OPT_A', "First Option", "Description one"),
            ('OPT_B', "Second Option", "Description two"),
        ),
        default='OPT_A',
    )

    def execute(self, context):
        return import_xcp(context, self.filepath, self.use_setting)


# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(ImportXCP.bl_idname, text="XCP Importer")


def register():
    bpy.utils.register_class(ImportXCP)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportXCP)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()

    # test call
    bpy.ops.importx.xcp('INVOKE_DEFAULT')