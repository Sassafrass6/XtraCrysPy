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