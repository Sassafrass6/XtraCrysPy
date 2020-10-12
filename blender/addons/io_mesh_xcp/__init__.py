bl_info = {
    "name": "XML Reader for XtraCrysPy",
    "blender": (2, 83, 0),
    "category": "Import-Export",
}

import bpy
from . import xcp_import
from . import xcp_panel

# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(xcp_import.ImportXCP.bl_idname, text="XCP Importer")


def register():
    bpy.utils.register_class(xcp_import.ImportXCP)
    bpy.utils.register_class(xcp_panel.EditorPanel)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(xcp_import.ImportXCP)
    bpy.utils.unregister_class(xcp_panel.EditorPanel)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()