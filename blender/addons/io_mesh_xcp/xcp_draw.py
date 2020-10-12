import bpy
import numpy as np

from . import xcp_math

def draw_ball(position, scale, name, collection, mtype="UV"):
    """draw a ball in blender

    Args:
        position (tuple(float, 3)): the position of the ball
        scale (float): radial scaling factor
        name (str): label to be displayed in blender outline
        collection (str): key for the collection to place new object in
        mtype (str, optional): method of drawing object. Defaults to "UV".

    Returns:
        bpy.types.Object: datablock for the created object
    """
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
    """draws a cylinder defined by two points

    Args:
        A (np.ndarray): 3-vector that specifies the location of one end face centre
        B (np.ndarray): 3-vector that specifies the locaiton of the other end face centre
        scale (float): the radial scaling for the cylinder
        name (str): label to be displayed in blender outline
        collection (str): key for the collection to place new object in
        mtype (str, optional): method of drawing object. Defaults to 'MESH'.

    Returns:
        bpy.types.Object: datablock for the created object
    """
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
    elif mtype == 'PYTHON':
        # mesh defined explicitly, taken from atomic plugin
        # should be faster than atomic plugin due to numpy and foreach_set (might want to verify this)
        sectors = 12 # tweak/move as needed
        length = np.linalg.norm(A - B)

        vertices = np.zeros((sectors * 2, 3))
        x = scale * np.cos(np.linspace(0, 2.0 * np.pi, sectors))
        y = scale * np.sin(np.linspace(0, 2.0 * np.pi, sectors))
        vertices[:sectors, 2] = 0.5 * length
        vertices[sectors:, 2] = -0.5 * length
        vertices[1:sectors, 0] = x[:-1]
        vertices[1:sectors, 1] = y[:-1]
        vertices[sectors+1:, 0] = x[:-1]
        vertices[sectors+1:, 1] = y[:-1]

        faces_side = []
        for i in range(sectors - 1):
            if i == sectors - 2:
                faces_side.append([i+1, 1, 1+sectors, i+1+sectors])
            else:
                faces_side.append([i+1, i+2, i+2+sectors, i+1+sectors])
        
        faces_caps = []
        for i in range(sectors - 1):
            if i == sectors - 2:
                face_top = [0, sectors-1, 1]
                face_bottom = [sectors, 2*sectors-1, sectors+1]
            else:
                face_top = [0, i+1, i+2]
                face_bottom = [sectors, sectors+1, sectors+2]
            faces_caps.append(face_top)
            faces_caps.append(face_bottom)
        
        

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

    if mtype == "MESH" or mtype == "NURBS" or mtype == "PYTHON":
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
        radiusA = bond["A"]["spinfo"]["scale"] / norm
        radiusB = bond["B"]["spinfo"]["scale"] / norm
        midpoint = (pointA * (1.0 + radiusB - radiusA) + pointB * (1.0 + radiusA - radiusB)) / 2.0
        labelA = bond["A"]["spinfo"]["label"]
        if label == labelA:
            pointB = midpoint
            species = bond["A"]["spinfo"]
        else:
            pointA = midpoint
            species = bond["B"]["spinfo"]
        n = (pointA - pointB) / norm
        gamma = np.dot(pointA, n)
        b = pointA - gamma * n
        n_b = b / np.linalg.norm(b)
        
        position = (pointA + pointB) / 2. - origin
        v = np.cross(n_b, n)
        v /= np.linalg.norm(v)
        p1 = position + n_b * 0.1
        p2 = position - n_b * 0.1
        p3 = position - v * 0.1
        p4 = position + v * 0.1

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
    stickobj = draw_stick(v1, v2, 0.2, "Stick.{}".format(label), collection_bonds_name, mtype='PYTHON')
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
