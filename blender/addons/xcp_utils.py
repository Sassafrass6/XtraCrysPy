import bpy

# copied from an older project and may not be up to date...
def removeAll(type=None):
    # possible type: 'MESH', 'CURVE', 'SURFACE', 'META', 'FONT', 'ARMATURE', 'LATTICE', 'EMPTY', 'CAMERA', 'LAMP'
    if type:
        bpy.ops.object.select_all(action='DESELECT')
        bpy.ops.object.select_by_type(type=type)
        bpy.ops.object.delete()
    # broken in blender 2.80
    else:
        override = bpy.context.copy()
        override['selected_objects'] = list(bpy.context.scene.objects)
        bpy.ops.object.delete(override)

def create_point(location):
    bpy.ops.object.empty_add(type="PLAIN_AXES", location=location)
    obj = bpy.context.active_object
    return obj

def track_to_constraint(obj, target):
    constraint = obj.constraints.new('TRACK_TO')
    constraint.target = target
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'

    return constraint

def limit_distance(obj, radius, target):
    constraint = obj.constraints.new('LIMIT_DISTANCE')
    constraint.target = target
    constraint.distance = radius
    constraint.limit_mode = "LIMITDIST_ONSURFACE"

def camera(origin=(0., 0., 0.), target=None, lens=45, clip_start=0.1, clip_end=200, type='PERSP', ortho_scale=6):
    bpy.ops.object.camera_add()
    camera = bpy.context.active_object
    camera.data.lens = lens
    camera.data.clip_start = clip_start
    camera.data.clip_end = clip_end
    camera.data.type = type
    camera.location = origin
    if type == 'ORTHO':
        camera.data.ortho_scale = ortho_scale

    if target: 
        track_to_constraint(camera, target)
    return camera

def lamp(origin=(0., 0., 0.), type='POINT', energy=1, color=(1,1,1), target=None):
    print('createLamp called')
    bpy.ops.object.light_add(type=type, location=origin)
    obj = bpy.context.active_object
    obj.data.energy = energy
    obj.data.color = color

    if target: 
        track_to_constraint(obj, target)
    return obj