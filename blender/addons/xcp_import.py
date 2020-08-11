bl_info = {
    "name": "XML Reader for XtraCrysPy",
    "blender": (2, 83, 0),
    "category": "Import-Export",
}

"""A WIP plugin for blender for reading xcp objects (XML-like)
Current design uses dictionaries over objects for now, since this is
a more direct interpretation of XML
"""
import bpy
import numpy as np # I think this comes with blender by default, doublecheck this
import xcp_io
import xcp_utils
import xcp_math
import xcp_draw

xcp_imports = {
    "xcp_io": xcp_io,
    "xcp_utils": xcp_utils,
    "xcp_math": xcp_math,
    "xcp_draw": xcp_draw
}

# should fix the import dependencies problem but currently doesn't work as intended
# When bpy is already in local, we know this is not the initial import...
if "bpy" in locals():
    # ...so we need to reload our submodule(s) using importlib
    import importlib
    for key, val in xcp_imports.items():
        if key in locals():
            importlib.reload(val)

def init_materials(species):
    for key in species:
        spec = species[key]
        spec["material"] = xcp_utils.create_material(key, spec["color"])

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
    surface = xcp_draw.draw_surface(spos, (cpos - mpos), radius * 3.0)
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
            xcp_draw.draw_duplivert_atom(xcp, atom_key_list, key, coll_atoms, def_origin, mtype=atom_draw_mode)
    else:
        for key, atom in xcp["ATOMS"].items():
            xcp["ATOMS"][key]["obj"] = xcp_draw.draw_atom(atom, key, coll_atoms, mtype=atom_draw_mode)
    # move into bond draw once we have methods to draw them in bulk
    if bond_dup:
        bond_map = xcp_utils.get_bond_map(xcp)
        for key, bond_key_list in bond_map.items():
            xcp_draw.draw_duplivert_bond(xcp, bond_key_list, key, coll_atoms, def_origin, default_materials["BOND"], bond_draw_mode, taper)
    else:
        for key in xcp["BONDS"]:
            bond = xcp["BONDS"][key]
            xcp["BONDS"][key]["obj"] = xcp_draw.draw_bond(bond, key, coll_atoms, default_materials["BOND"], bond_draw_mode, taper)
    if "CAMERA" in xcp["SCENE"]:
        # untested
        camera = xcp["SCENE"]["CAMERA"]
        cameraobj = xcp_utils.camera_from_input(camera["position"], camera["rotation"], coll_atoms)
    if "FRAME" in xcp:
        frameobj = xcp_draw.draw_frame(xcp["FRAME"], coll_frame)
        
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