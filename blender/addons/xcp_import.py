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
import xcp_math

def add_camera(position, rotation, collection):
    cam = bpy.data.cameras.new("Camera")
    cam.lens = 18

    # create the first camera object
    objref = bpy.data.objects.new("Camera", cam)
    objref.location = position
    objref.rotation_euler = rotation
    bpy.data.collections[collection].objects.link(objref)
    return objref

def draw_ball(position, scale, name, collection, mtype="UV"):
    # TODO add this to a collection of sorts
    if mtype == "UV":
        bpy.ops.mesh.primitive_uv_sphere_add(
            segments=64, ring_count=48, 
            align='WORLD', enter_editmode=False)
    elif mtype == "META":
        bpy.ops.object.metaball_add(type='BALL', align='WORLD')
    elif mtype == "NURBS":
        bpy.ops.surface.primitive_nurbs_surface_sphere_add(
            align='WORLD'
        )

    ball = bpy.context.view_layer.objects.active     
    bpy.ops.collection.objects_remove_all()
    bpy.data.collections[collection].objects.link(ball)
    ball.location = position
    ball.scale = (scale, scale, scale)
    ball.name = name

    return ball

def draw_stick(A, B, scale, name, collection, mtype='MESH'):
    # next, get rotation matrix
    position = (A + B) / 2.0
    rotator = xcp_math.vec_to_rotator(A - B)    
    quat = xcp_math.rotator_to_quaternion(rotator)

    # call create mesh
    if mtype == 'MESH':
        bpy.ops.mesh.primitive_cylinder_add(align='WORLD')
    elif mtype == 'META':
        bpy.ops.object.metaball_add(type='CAPSULE', align='WORLD')
        meta_mod = 1 / (np.sqrt(2))
        quat = xcp_math.quaternion_multiply((meta_mod, 0., meta_mod, 0.), quat)
    elif mtype == 'NURBS':
        bpy.ops.surface.primitive_nurbs_surface_cylinder_add(align='WORLD')
    else:
        return

    stick = bpy.context.view_layer.objects.active
    bpy.ops.collection.objects_remove_all()
    bpy.data.collections[collection].objects.link(stick)
    stick.location = tuple(position)
    stick.rotation_mode = "QUATERNION"
    stick.rotation_quaternion = quat

    norm = np.linalg.norm(A - B)
    if mtype == 'MESH' or mtype == 'NURBS':
        stick.scale = (scale, scale, norm / 2.0)
    elif mtype == 'META':
        stick.scale = (norm / 4.0, scale, scale)
    stick.name = name

    return stick

def draw_bezier(A, B, scale, name, collection, taper):
    midpoint = (A + B) / 2.0
    bpy.ops.curve.primitive_bezier_curve_add(radius=1.0,
        location=(0.0, 0.0, 0.0))

    curve = bpy.context.view_layer.objects.active
    curve.location = midpoint
    bezier_points = curve.data.splines[0].bezier_points
    bezier_points[0].co = A - midpoint
    bezier_points[0].handle_left = (A - midpoint) * 1.1
    bezier_points[0].handle_right = (0., 0., 0.)
    bezier_points[1].co = B - midpoint
    bezier_points[1].handle_left = (0., 0., 0.)
    bezier_points[1].handle_right = (B - midpoint) * 1.1
    
    curve.data.bevel_depth = scale
    curve.data.bevel_resolution = 10
    if taper:
        curve.data.taper_object = taper
    bpy.ops.object.modifier_add(type="SOLIDIFY")
    curve.modifiers["Solidify"].thickness = 0.1

    return curve

def draw_surface(center, norm, scale):
    rotator = xcp_math.vec_to_rotator(norm)
    quat = xcp_math.rotator_to_quaternion(rotator)
    bpy.ops.mesh.primitive_plane_add(location=center)
    surface = bpy.context.view_layer.objects.active
    surface.rotation_mode = "QUATERNION"
    surface.rotation_quaternion = quat
    surface.scale = (scale, scale, scale)
    return surface

def draw_atom(atom, uid, collection, mtype='UV'):
    if mtype == "META":
        name = ""
    else:
        name = "Atom.{}.{:03d}.{}".format(collection[-3:], uid, atom["spinfo"]["label"])
    objref = draw_ball(atom["position"], atom["spinfo"]["scale"], name, collection, mtype)
    objref.data.materials.append(atom["spinfo"]["material"])
    return objref

def draw_duplivert_atom(data, key_list, label, collection, origin, mtype='UV'):
    # mostly taken from the official pdb addon for drawing duplivert spheres
    vertices = []
    species = None
    for key in key_list:
        atom = data["ATOMS"][key]
        # all of these should be the same, should also sanity check that this is not None
        species = atom["spinfo"]
        vertices.append(atom["position"] - origin)
    
    if species is None:
        print("Species info not linked to atom in xml input!")
    
    collection_species_name = "AtomicSpecies.{}".format(label)
    collection_species = bpy.data.collections.new(collection_species_name)
    bpy.data.collections[collection].children.link(collection_species)

    # this mesh contains the information about all the duplicates
    atom_dup_mesh = bpy.data.meshes.new("AtomDupMesh.{}".format(label))
    atom_dup_mesh.from_pydata(vertices, [], [])
    atom_dup_mesh.update()
    # this is the main mesh that, when linked, will correctly display the duplicates
    atom_mesh = bpy.data.objects.new("AtomMesh.{}".format(label), atom_dup_mesh)
    collection_species.objects.link(atom_mesh)

    ballobj = draw_ball((0., 0., 0.), species["scale"], "Ball.{}".format(label), collection_species_name, mtype)
    # if not metaball, hide the original!
    if mtype != "META":
        ballobj.hide_set(True)
    ballobj.data.materials.append(species["material"])
    ballobj.parent = atom_mesh

    atom_mesh.instance_type = "VERTS"
    atom_mesh.location = origin

    return atom_mesh, collection_species

def draw_bond(bond, uid, collection, default_material, mtype="MESH", taper=None):
    if mtype == "META":
        name = ""
    else:
        name = "Bond.{}.{:03d}".format(collection[-3:], uid)
    # in the case for multiple bond objects
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

    if mtype == "MESH" or mtype == "NURBS":
        objref.append(draw_stick(pointA, midpoint, 0.2, nameA, collection, mtype))
        objref.append(draw_stick(midpoint, pointB, 0.2, nameB, collection, mtype))
        objref[0].data.materials.append(bond["A"]["spinfo"]["material"])
        objref[1].data.materials.append(bond["B"]["spinfo"]["material"])
    elif mtype == "META":
        objref.append(draw_stick(pointA, pointB, 0.2, name, collection, mtype='META'))
        objref[0].data.materials.append(default_material)
    elif mtype == "BEZIER":
        vecA = (pointB - pointA) / norm
        vecB = (pointA - pointB) / norm
        pointA += bond["A"]["spinfo"]["scale"] * vecA * 0.9
        pointB += bond["B"]["spinfo"]["scale"] * vecB * 0.9
        objref.append(draw_bezier(pointA, pointB, 0.4, name, collection, taper))
        objref[0].data.materials.append(default_material)
    return objref

def draw_duplivert_bond(data, key_list, label, collection, origin, default_material, mtype="MESH", taper=None):
    # mostly taken from the official pdb addon for drawing duplivert sticks
    vertices = []
    faces = []
    i = 0
    species = None
    for key in key_list:
        bond = data["BONDS"][key]   
        pointA = np.array(bond["A"]["position"])
        pointB = np.array(bond["B"]["position"])
        norm = np.linalg.norm(pointA - pointB)
        n = (pointA - pointB) / norm
        radiusA = bond["A"]["spinfo"]["scale"] / norm
        radiusB = bond["B"]["spinfo"]["scale"] / norm
        midpoint = (pointA * (1.0 + radiusB - radiusA) + pointB * (1.0 + radiusA - radiusB)) / 2.0
        labelA = bond["A"]["spinfo"]["label"]
        if label == labelA:
            point = pointA
            species = bond["A"]["spinfo"]
        else:
            point = pointB
            species = bond["B"]["spinfo"]
        
        position = (point + midpoint) / 2. - origin
        if n[0] > 0.9:
            v1 = np.array([0., 1., 0.])
        else:
            v1 = np.array([1., 0., 0.])
        v2 = np.cross(n, v1)
        v3 = np.cross(n, v2)
        p1 = position + v2 / np.linalg.norm(v2) * 0.2
        p2 = position - v2 / np.linalg.norm(v2) * 0.2
        p3 = position + v3 / np.linalg.norm(v3) * 0.2
        p4 = position - v3 / np.linalg.norm(v3) * 0.2

        vertices.append(p1)
        vertices.append(p2)
        vertices.append(p3)
        vertices.append(p4)

        faces.append((i*4+0, i*4+2, i*4+1, i*4+3))
        i += 1

    collection_bonds_name = "Bonds.{}".format(label)
    collection_bonds = bpy.data.collections.new(collection_bonds_name)
    bpy.data.collections[collection].children.link(collection_bonds)

    # this mesh contains the information about all the duplicates
    bond_dup_mesh = bpy.data.meshes.new("BondDupMesh.{}".format(label))
    bond_dup_mesh.from_pydata(vertices, [], faces)
    bond_dup_mesh.update()
    # this is the main mesh that, when linked, will correctly display the duplicates
    bond_mesh = bpy.data.objects.new("BondMesh.{}".format(label), bond_dup_mesh)
    collection_bonds.objects.link(bond_mesh)

    v1 = np.array([0., 0., 0.])
    v2 = np.array([0., 1., 0.])
    stickobj = draw_stick(v1, v2, 0.2, "Stick.{}".format(label), collection_bonds_name, mtype)
    stickobj.hide_set(True)
    stickobj.data.materials.append(species["material"])
    stickobj.parent = bond_mesh

    bond_mesh.instance_type = "FACES"
    bond_mesh.location = origin

    return bond_mesh, collection_bonds


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
    cameraobj = xcp_utils.camera(target=constraint)
    xcp_utils.limit_distance(cameraobj, radius, constraint)
    l1pos = origin + np.array((1., 1., 0.))
    l2pos = origin + np.array((-1., -1., -1.))
    l3pos = origin + np.array((0., 1., 1.))
    lampobj1 = xcp_utils.lamp(origin=l1pos, energy=10 * radius ** 2, target=constraint)
    lampobj2 = xcp_utils.lamp(origin=l2pos, energy=10 * radius ** 2, target=constraint)
    lampobj3 = xcp_utils.lamp(origin=l3pos, energy=10 * radius ** 2, target=constraint)
    xcp_utils.limit_distance(lampobj1, radius - 4.0, constraint)
    xcp_utils.limit_distance(lampobj2, radius - 4.0, constraint)
    xcp_utils.limit_distance(lampobj3, radius - 4.0, constraint)

    return (constraint, [cameraobj], [lampobj1, lampobj2, lampobj3])

def add_default_background(camera_ref, molecule_ref, radius):
    cpos = np.array(camera_ref.location)
    mpos = np.array(molecule_ref.location)
    spos = cpos + 2.5 * (mpos - cpos)
    surface = draw_surface(spos, (cpos - mpos), radius * 3.0)
    return surface

def add_sphere_background(origin, radius):
    bpy.ops.surface.primitive_nurbs_surface_sphere_add()
    surface = bpy.context.view_layer.objects.active 
    surface.location = origin
    surface.scale = (radius, radius, radius)
    surface.hide_viewport = True
    return surface

def import_xcp(context, filepath, clear_world, default_view, join_mode, atom_draw_mode, bond_draw_mode,
        atom_dup, bond_dup):
    if clear_world:
        xcp_utils.removeAll()
    xcp = xcp_io.read_xcp(filepath)
    def_cameras = []
    def_lights = []
    def_origin = get_center(xcp["ATOMS"])

    if default_view:
        def_radius = get_span(xcp["ATOMS"]) * 2.0
        point, def_cameras, def_lights = add_default_view(def_origin, def_radius)
        #def_background = add_default_background(def_cameras[0], point, def_radius)
        def_background = add_sphere_background(def_origin, def_radius * 1.5)

    coll_id = "{:03d}".format(get_collection_id() + 1)
    coll_atoms = "AtomGroup.{}".format(coll_id)
    coll_frame = "CellFrame.{}".format(coll_id)
    collref_atoms = bpy.data.collections.new(coll_atoms) 
    collref_frame = bpy.data.collections.new(coll_frame) 
    bpy.context.scene.collection.children.link(collref_atoms)
    bpy.context.scene.collection.children.link(collref_frame)

    # do some processing in between perhaps and then finally here, draw every atom
    init_materials(xcp["SPECIES"])
    default_materials = xcp_utils.init_default_materials()
    if bond_draw_mode == "BEZIER":
        taper = xcp_utils.create_taper(label="Taper")
    else:
        taper = None

    if atom_dup:
        species_map = xcp_utils.get_species_map(xcp)
        for key, atom_key_list in species_map.items():
            draw_duplivert_atom(xcp, atom_key_list, key, coll_atoms, def_origin, mtype=atom_draw_mode)
    else:
        for key, atom in xcp["ATOMS"].items():
            xcp["ATOMS"][key]["obj"] = draw_atom(atom, key, coll_atoms, mtype=atom_draw_mode)
    # move into bond draw once we have methods to draw them in bulk
    if bond_dup:
        bond_map = xcp_utils.get_bond_map(xcp)
        for key, bond_key_list in bond_map.items():
            draw_duplivert_bond(xcp, bond_key_list, key, coll_atoms, def_origin, default_materials["BOND"], bond_draw_mode, taper)
    else:
        for key in xcp["BONDS"]:
            bond = xcp["BONDS"][key]
            xcp["BONDS"][key]["obj"] = draw_bond(bond, key, coll_atoms, default_materials["BOND"], bond_draw_mode, taper)
    if "CAMERA" in xcp["SCENE"]:
        # untested
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
            ('OPT_C', "Molecule Join", "Entire molecule is joined"),
        ),
        default='OPT_A',
    )

    atom_draw_mode: EnumProperty(
        name="Atom Representation",
        description="Determines how the atoms are displayed.",
        items=(
            ('UV', "Mesh", "Standard molecular drawing"),
            ('META', "Metaballs", "Each object represented by a metaball"),
            ('NURBS', "NURBS", "Better for rendering, but more expensive")
        )
    )

    atom_duplivert: BoolProperty(
        name="Duplivert (single atomic species)",
        description="Draws all atoms of same species as duplivert",
        default=False,
    )

    bond_draw_mode: EnumProperty(
        name="Bond Representation",
        description="Determines how the bonds are displayed.",
        items=(
            ('MESH', "Mesh", "Standard molecular drawing"),
            ('META', "Metaballs", "Each object represented by a metaball"),
            ('BEZIER', "Plastic", "Bonds are represented as bezier curves"),
            ('NURBS', "NURBS", "Better for rendering, but more expensive")
        )
    )

    bond_duplivert: BoolProperty(
        name="Duplivert (bonds)",
        description="Draws all bonds of same species as duplivert",
        default=False,
    )

    def execute(self, context):
        return import_xcp(context, self.filepath, self.clear_world, 
            self.default_view, self.join_mode, self.atom_draw_mode, self.bond_draw_mode,
            self.atom_duplivert, self.bond_duplivert)


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