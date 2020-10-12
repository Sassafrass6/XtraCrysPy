from bpy.types import Panel

class EditorPanel(Panel):
    bl_label = "XtraCrysPy Blender Utilities"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Create"
    bl_idname = "XCPUtils"

    # Guarantee that panel only opened when addon is activated
    @classmethod
    def poll(cls, context):
        return context.preferences.addons[__package__].preferences.bool_utility

    def draw(self, context):
        layout = self.layout
        