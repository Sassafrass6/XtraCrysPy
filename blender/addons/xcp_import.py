bl_info = {
    "name": "XML Reader for XtraCrysPy",
    "blender": (2, 83, 0),
    "category": "Import-Export",
}

"""A WIP plugin for blender for reading xcp objects (XML-like)
Current design uses dictionaries over objects for now, since this is
a more direct interpretation of XML
"""

# should fix the import dependencies problem but currently doesn't work as intended
# When bpy is already in local, we know this is not the initial import...
if "bpy" in locals():
    # ...so we need to reload our submodule(s) using importlib
    import importlib
    if "xcp_io" in locals():
        importlib.reload(xcp_io)
    if "xcp_utils" in locals():
        importlib.reload(xcp_utils)

import bpy
import numpy as np # I think this comes with blender by default, doublecheck this
import xcp_io
import xcp_utils

def add_camera(position, rotation, collection):
    cam = bpy.data.cameras.new("Camera")
    cam.lens = 18

    # create the first camera object
    objref = bpy.data.objects.new("Camera", cam)
    objref.location = position
    objref.rotation_euler = rotation
    bpy.data.collections[collection].objects.link(objref)
    return objref

def draw_ball(position, scale, name, collection):
    # TODO add this to a collection of sorts
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=64, ring_count=48, 
        align='WORLD', enter_editmode=False, 
        location=(0., 0., 0.), rotation=(0., 0., 0.))
        
    ball = bpy.context.view_layer.objects.active
    bpy.ops.collection.objects_remove_all()
    bpy.data.collections[collection].objects.link(ball)
    ball.location = position
    ball.scale = (scale, scale, scale)
    ball.name = name

    return ball

def draw_stick(A, B, scale, name, collection):
    # next, get rotation matrix
    origin = (A + B) / 2.0
    vec1 = A - B
    norm = np.linalg.norm(vec1)
    vec1 /= norm
    vec2 = np.array([0, 0, 1])
    if (vec1 == vec2).all():
        rotator = np.identity(3)
    elif (vec1 == -vec2).all():
        rotator = np.identity(3) * -1.0
        rotator[0, 0] += 2.0
    else:
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
    bpy.ops.collection.objects_remove_all()
    bpy.data.collections[collection].objects.link(stick)

    # scale accordingly
    scale4 = np.eye(4, dtype=float)
    scale4[0, 0] = scale
    scale4[1, 1] = scale
    scale4[2, 2] = norm / 2

    # apply these transformations
    stick.data.transform(scale4)
    stick.data.transform(rot4)
    stick.data.update()
    stick.name = name

    return stick

def draw_atom(atom, uid, collection):
    name = "Atom.{}.{:03d}.{}".format(collection[-3:], uid, atom["spinfo"]["label"])
    objref = draw_ball(atom["position"], atom["spinfo"]["scale"], name, collection)
    return objref

def draw_bond(bond, uid, collection):
    nameA = "Bond.{}.{:03d}.A".format(collection[-3:], uid)
    nameB = "Bond.{}.{:03d}.B".format(collection[-3:], uid)
    pointA = np.array(bond["A"]["position"])
    pointB = np.array(bond["B"]["position"])
    norm = np.linalg.norm(pointA - pointB)
    radiusA = bond["A"]["spinfo"]["scale"] / norm
    radiusB = bond["B"]["spinfo"]["scale"] / norm
    # this equation takes in account the radii of each atom and sets the midpoint to reveal
    # an equal portion of bond
    midpoint = (pointA * (1.0 + radiusB - radiusA) + pointB * (1.0 + radiusA - radiusB)) / 2.0
    objref = []

    objref.append(draw_stick(pointA, midpoint, 0.2, nameA, collection))
    objref.append(draw_stick(midpoint, pointB, 0.2, nameB, collection))
    return objref

def draw_frame(frame, collection):
    # TODO move defaults into single location
    frame_scale = 0.1
    objref = []
    for key, vertex in frame["VERTICES"].items():
        name = "FrameVertex.{:03d}".format(key)
        objref.append(draw_ball(vertex["position"], frame_scale, name, collection))
    for key, edge in frame["EDGES"].items():
        name = "FrameEdge.{:03d}".format(key)
        pointA = np.array(edge["A"]["position"])
        pointB = np.array(edge["B"]["position"])
        objref.append(draw_stick(pointA, pointB, frame_scale, name, collection))
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

def get_collection_id():
    id = -1
    sub_collections = bpy.context.scene.collection.children.keys()
    for sc in sub_collections:
        if sc.startswith("AtomGroup"):
            id = max(id, int(sc.split(".")[1]))
    return id

def get_center(atoms):
    center = np.zeros(3)
    for atom in atoms.values():
        center += np.array(atom["position"])
    center /= len(atoms)
    return center

def get_span(atoms):
    pmin = np.zeros(3)
    pmax = np.zeros(3)
    for atom in atoms.values():
        pmin = np.minimum(pmin, np.array(atom["position"]))
        pmax = np.maximum(pmax, np.array(atom["position"]))
    return np.max(pmax - pmin)

def add_default_view(origin, radius):
    constraint = xcp_utils.create_point(origin)
    # temporary
    cameraobj = xcp_utils.camera(target=constraint)
    xcp_utils.limit_distance(cameraobj, radius, constraint)
    l1pos = origin + np.array((1., 1., 0.))
    l2pos = origin + np.array((-1., -1., -1.))
    l3pos = origin + np.array((0., 1., 1.))
    lampobj1 = xcp_utils.lamp(origin=l1pos, energy=10000, target=constraint)
    lampobj2 = xcp_utils.lamp(origin=l2pos, energy=10000, target=constraint)
    lampobj3 = xcp_utils.lamp(origin=l3pos, energy=10000, target=constraint)
    xcp_utils.limit_distance(lampobj1, radius - 4.0, constraint)
    xcp_utils.limit_distance(lampobj2, radius - 4.0, constraint)
    xcp_utils.limit_distance(lampobj3, radius - 4.0, constraint)

def import_xcp(context, filepath, clear_world, default_view, join_mode):
    if clear_world:
        xcp_utils.removeAll()
    xcp = xcp_io.read_xcp(filepath)

    if default_view:
        def_origin = get_center(xcp["ATOMS"])
        def_radius = get_span(xcp["ATOMS"]) * 2.0
        add_default_view(def_origin, def_radius)

    coll_id = "{:03d}".format(get_collection_id() + 1)
    coll_atoms = "AtomGroup.{}".format(coll_id)
    coll_frame = "CellFrame.{}".format(coll_id)
    collref_atoms = bpy.data.collections.new(coll_atoms) 
    collref_frame = bpy.data.collections.new(coll_frame) 
    bpy.context.scene.collection.children.link(collref_atoms)
    bpy.context.scene.collection.children.link(collref_frame)

    # do some processing in between perhaps and then finally here, draw every atom
    init_materials(xcp["SPECIES"])
    for key in xcp["ATOMS"]:
        atom = xcp["ATOMS"][key]
        atomobj = draw_atom(atom, key, coll_atoms)
        # link materials (expects init_materials to have been called)
        atomobj.data.materials.append(atom["spinfo"]["material"])
    for key in xcp["BONDS"]:
        bond = xcp["BONDS"][key]
        bondobj = draw_bond(bond, key, coll_atoms)
        bondobj[0].data.materials.append(bond["A"]["spinfo"]["material"])
        bondobj[1].data.materials.append(bond["B"]["spinfo"]["material"])
    if "CAMERA" in xcp["SCENE"]:
        camera = xcp["SCENE"]["CAMERA"]
        cameraobj = add_camera(camera["position"], camera["rotation"], coll_atoms)
    if "FRAME" in xcp:
        frameobj = draw_frame(xcp["FRAME"], coll_frame)
        
    # tells blender that the addon main function is successful
    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator


class ImportXCP(Operator, ImportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "importx.xml"
    bl_label = "Import XML Data Generated by XtraCrysPy"

    # ImportHelper mixin class uses this
    filename_ext = ".xml"

    filter_glob: StringProperty(
        default="*.xml",
        options={'HIDDEN'},
        maxlen=255,  # Max internal buffer length, longer would be clamped.
    )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    clear_world: BoolProperty(
        name="Clear World",
        description="Removes everything from current world",
        default=False,
    )

    default_view: BoolProperty(
        name="Default View",
        description="Adds a camera and lights defaulted to look at the molecule at some sensible distance",
        default=False,
    )

    # TODO join mode should remove inside faces but I'm not certain how to do this
    join_mode: EnumProperty(
        name="Join Mode",
        description="Determines whether meshes are joined. Joining facilitates control.",
        items=(
            ('OPT_A', "No Join", "Each atom and bond is a separate mesh"),
            ('OPT_B', "Atom Join", "Bonds are joined to each respective atom"),
            ('OPT_C', "Molecule Join", "Entire molecule is joined")
        ),
        default='OPT_A',
    )

    def execute(self, context):
        return import_xcp(context, self.filepath, self.clear_world, self.default_view, self.join_mode)


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