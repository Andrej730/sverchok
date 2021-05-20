# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####
import bpy
from bpy.props import EnumProperty, FloatProperty, FloatVectorProperty, BoolProperty
import numpy as np
from sverchok.node_tree import SverchCustomTreeNode
from sverchok.data_structure import updateNode, list_match_func, numpy_list_match_modes, numpy_list_match_func
from sverchok.utils.sv_itertools import recurse_f_level_control


def scale_meshes(params, constant, matching_f):
    '''
    params are verts, movement vectors, strength
    - verts and  movement vectors should be list as [[[float, float, float],],] (Level 3)
    - Strength should be list as [[float, float, ..], [float, ..], ..] (Level 2)
    desired_levels = [3, 3, 2]
    constant are the function options (data that does not need to be matched)
    matching_f stands for list matching formula to use
    '''
    result = []
    match_mode, output_numpy = constant
    params = matching_f(params)

    local_match = numpy_list_match_func[match_mode]

    for props in zip(*params):
        verts, centers, scale, strength = local_match([np.array(p) for p in props])
        verts_out = centers + (verts - centers) * scale * strength[:, np.newaxis]
        result.append(verts_out if output_numpy else verts_out.tolist())

    return result


class SvScaleNodeMk3(bpy.types.Node, SverchCustomTreeNode):
    """
    Triggers: Scale vertices
    Tooltip: Scales vectors from a center point

    """

    bl_idname = 'SvScaleNodeMk3'
    bl_label = 'Scale'
    bl_icon = 'ORIENTATION_VIEW'
    sv_icon = 'SV_SCALE'


    centers: FloatVectorProperty(
        name='Centers', description='Center of the scaling transform',
        size=3, default=(0, 0, 0),
        update=updateNode)
    scale: FloatVectorProperty(
        name='Scale', description='Center of the scaling transform',
        size=3, default=(1, 1, 1),
        update=updateNode)

    multiplier: FloatProperty(
        name='Multiplier', description='Multiplier factor',
        default=1.0, update=updateNode)

    list_match: EnumProperty(
        name="List Match",
        description="Behavior on different list lengths",
        items=numpy_list_match_modes, default="REPEAT",
        update=updateNode)

    output_numpy: BoolProperty(
        name='Output NumPy',
        description='Output NumPy arrays',
        default=False, update=updateNode)

    def sv_init(self, context):

        self.inputs.new('SvVerticesSocket', 'Vertices')
        self.inputs.new('SvVerticesSocket', 'Centers').prop_name = 'centers'
        self.inputs.new('SvVerticesSocket', 'Scale').prop_name = 'scale'
        self.inputs.new('SvStringsSocket', 'Strength').prop_name = 'multiplier'


        self.outputs.new('SvVerticesSocket', 'Vertices')

    def migrate_from(self, old_node):
        self.multiplier = old_node.factor_

    def draw_buttons_ext(self, context, layout):
        '''draw buttons on the N-panel'''
        layout.prop(self, 'output_numpy')
        layout.prop(self, 'list_match', expand=False)

    def rclick_menu(self, context, layout):
        layout.prop(self, 'output_numpy')
        layout.prop_menu_enum(self, "list_match", text="List Match")

    def process(self):
        inputs, outputs = self.inputs, self.outputs

        if not outputs[0].is_linked:
            return

        result = []

        params = [si.sv_get(default=[[]], deepcopy=False) for si in inputs]

        matching_f = list_match_func[self.list_match]
        desired_levels = [3, 3, 3, 2]
        ops = [self.list_match, self.output_numpy]

        result = recurse_f_level_control(params, ops, scale_meshes, matching_f, desired_levels)

        self.outputs[0].sv_set(result)

classes = [SvScaleNodeMk3]
register, unregister = bpy.utils.register_classes_factory(classes)
