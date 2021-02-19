bl_info = {
    "name": "XML Reader for XtraCrysPy",
    "blender": (2, 83, 0),
    "category": "Import-Export",
}

import bpy
from bpy.props import PointerProperty
from . import xcp_import
from . import xcp_panel

# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(xcp_import.ImportXCP.bl_idname, text="XCP Importer")

classes = (
    xcp_import.ImportXCP,
    xcp_panel.PanelProperties,
    xcp_panel.XCP_OT_Material,
    xcp_panel.XCP_PT_Material
)

def register():
    for cls in classes:
        print("add {}".format(cls))
        bpy.utils.register_class(cls)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)
    bpy.types.Scene.xcptool = PointerProperty(type=xcp_panel.PanelProperties)


def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)
    del bpy.types.Scene.xcptool


if __name__ == "__main__":
    register()