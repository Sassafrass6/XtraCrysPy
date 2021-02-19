from bpy.props import EnumProperty
from bpy.types import Operator, Panel, PropertyGroup
import bpy

from . import xcp_data
class XCP_PT_Material(Panel):
    bl_idname = "XCP_PT_material"
    bl_label = "XCP Material Panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "XCP"

    @classmethod
    def poll(self, context):
        return context.object is not None

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        tool = scene.xcptool
        layout.label(text="Material Editor")

        layout.prop(tool, "material_choice")
        layout.operator("xcp.material")
        layout.separator()


class XCP_OT_Material(Operator):
    bl_label = "Update"
    bl_idname = "xcp.material"

    def execute(self, context):
        # TODO consider making all the materials on startup and then cloning and setting specific colors etc. 

        scene = context.scene
        tool = scene.xcptool

        # print the values to the console
        select_list = bpy.context.selected_objects
        material_list = []
        for s in select_list:
            for m in s.data.materials:
                material_list.append(m)
        print(material_list)

        for m in material_list:
            color = None
            # check if the pre-existing material setup aligns with the expected 
            if not m.use_nodes:
                # if not, try to get the color and set `use_nodes` to true, to be reset later
                color = m.diffuse_color
                m.use_nodes = True

            nodes = m.node_tree.nodes
            links = m.node_tree.links

            # if default setup, get the color from the principled bsdf, otherwise (for now)
            # complain and continue with no color
            # TODO sort out colors in non-standard cases
            # TODO maybe the atom color needs to be set less ambigiously?
            if "Principled BSDF" in nodes:
                color = nodes["Principled BSDF"].inputs[0].default_value
            if color is None:
                print("Unable to figure out the initial color for {}, setting to a default for the user to change".format(m))
                color = (0., 0., 0.)

            # recursively clear the existing `node_tree`
            links.clear()
            nodes.clear()

            # read in the data and build the tree
            material_data = xcp_data.material_data[tool.material_choice]
            new_node_dict = {}
            for node_key in material_data:
                node_name = "ShaderNode{}".format(node_key.split('.')[0])
                node_data = material_data[node_key]
                node = nodes.new(node_name)
                new_node_dict[node_key] = node
                # set all the defaults
                if 'inputs' in node_data:
                    for i, data in enumerate(node_data['inputs']):
                        if data == None:
                            # do nothing
                            continue
                        if type(data) is xcp_data.Link:
                            # try to form a link (it's ok if this fails as the other connected node should
                            # try to connect back)
                            if data.node not in new_node_dict:
                                continue
                            links.new(node.inputs[i], new_node_dict[data.node].outputs[data.target_index])
                        # otherwise try to just set the value
                        else:
                            node.inputs[i].default_value = data
                if 'outputs' in node_data:
                    for i, data in enumerate(node_data['outputs']):
                        if data == None:
                            # do nothing
                            continue
                        if type(data) is xcp_data.Link:
                            # try to form a link (it's ok if this fails as the other connected node should
                            # try to connect back)
                            if data.node not in new_node_dict:
                                continue
                            links.new(new_node_dict[data.node].inputs[data.target_index], node.outputs[i])
                for key in node_data:
                    if key == 'inputs' or key == 'outputs':
                        continue
                    # set anything else
                    setattr(node, key, node_data[key])

        return {'FINISHED'}

class PanelProperties(PropertyGroup):
    material_choice: EnumProperty(
        name="Change Material",
        description="Select a predefined material for the object",
        items=[
            ('mat.glossy', "Glossy", ""),
            ('mat.opaque', "Opaque", ""),
            ('mat.metallic', "Metal", ""),
            ('mat.glass', "Glass", ""),
            ('mat.translucent', "Translucent", ""),
            ('mat.plastic1', "Plastic 1", ""),
            ('mat.plastic2', "Plastic 2", ""),
            ('mat.chalk', "Chalk", ""),
            # ('mat.sandstone', "Sandstone", ""),
            # ('mat.granite', "Granite", "")
        ]
    )
