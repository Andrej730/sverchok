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
from bpy.props import EnumProperty, FloatVectorProperty, IntProperty, FloatProperty, IntVectorProperty, BoolProperty, StringProperty, CollectionProperty, PointerProperty

from sverchok.node_tree import SverchCustomTreeNode
from sverchok.data_structure import updateNode, multi_socket, ensure_min_nesting, ensure_nesting_level
from sverchok.utils.nodes_mixins.show_3d_properties import Show3DProperties
from sverchok.utils.sv_logging import sv_logger
from sverchok.utils.handle_blender_data import correct_collection_length
from sverchok.utils.sv_operator_mixins import SvGenericNodeLocator
import numpy as np
import re
from itertools import compress, chain
from mathutils import Vector, Quaternion, Euler, Matrix
import json

def ShowMessageBox(message = "", title = "Message Box", icon = 'INFO'):
    def draw(self, context):
        self.layout.label(text=message)
    bpy.context.window_manager.popup_menu(draw, title = title, icon = icon)

class SvListInputBoolEntry(bpy.types.PropertyGroup):
    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: BoolProperty(
        name = "Boolean",
        default = True,
        update=update_entry,
    ) # type: ignore

class SvListInputIntEntry(bpy.types.PropertyGroup):
    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: IntProperty(
        name = "Integer number",
        default = 0,
        update=update_entry,
    ) # type: ignore

class SvListInputFloatEntry(bpy.types.PropertyGroup):

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value
        
    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: FloatProperty(
        name = "Float number",
        default = 0,
        #subtype='ANGLE',
        update=update_entry,
    ) # type: ignore
    
    NONE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='NONE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    PERCENTAGE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='PERCENTAGE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    FACTOR: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='FACTOR',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    ANGLE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='ANGLE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    TIME: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='TIME',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    TIME_ABSOLUTE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='TIME_ABSOLUTE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    DISTANCE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='DISTANCE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    DISTANCE_CAMERA: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='DISTANCE_CAMERA',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    TEMPERATURE: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='TEMPERATURE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore
    
    # WAVELENGTH: FloatProperty(
    #     name = "Float number",
    #     default = 0,
    #     subtype='WAVELENGTH',
    #     get = get_elem, set = set_elem,
    #     update=update_entry,
    # ) # type: ignore
    
    POWER: FloatProperty(
        name = "Float number",
        default = 0,
        subtype='POWER',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

class SvListInputVectorEntry(bpy.types.PropertyGroup):

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value

    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: FloatVectorProperty(
        name = "Vector",
        default = (0., 0., 0.,),
        #subtype='COLOR', 
        #subtype='XYZ',
        update=update_entry,
    ) # type: ignore

    NONE: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='NONE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    COLOR: FloatVectorProperty(
        name = "Color",
        default = (1., 1., 1.,),
        min = 0.0,
        max = 1.0,
        subtype='COLOR',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    TRANSLATION: FloatVectorProperty(
        name = "Translation",
        default = (0., 0., 0.,),
        subtype='TRANSLATION',
        get = get_elem, set = set_elem
    ) # type: ignore

    DIRECTION: FloatVectorProperty(
        name = "Float number",
        default = (1., 1., 1.,),
        subtype='DIRECTION',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    VELOCITY: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='VELOCITY',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    ACCELERATION: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='ACCELERATION',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    # MATRIX: FloatVectorProperty(
    #     name = "Float number",
    #     default = 0,
    #     subtype='MATRIX',
    #     get = get_elem, set = set_elem
    # ) # type: ignore

    EULER: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='EULER',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    # QUATERNION: FloatVectorProperty(
    #     name = "Float number",
    #     default = (0., 0., 0.,),
    #     subtype='QUATERNION',
    #     get = get_elem, set = set_elem,
    #     update=update_entry,
    # ) # type: ignore

    AXISANGLE: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='AXISANGLE',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    XYZ: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='XYZ',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    XYZ_LENGTH: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='XYZ_LENGTH',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    COLOR_GAMMA: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='COLOR_GAMMA',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    COORDINATES: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='COORDINATES',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    LAYER: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='LAYER',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    LAYER_MEMBER: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0.,),
        subtype='LAYER_MEMBER',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

quaternion_euler_orders = [
    ('XYZ', 'XYZ', "", 0),
    ('XZY', 'XZY', "", 1),
    ('YXZ', 'YXZ', "", 2),
    ('YZX', 'YZX', "", 3),
    ('ZXY', 'ZXY', "", 4),
    ('ZYX', 'ZYX', "", 5)
]

quaternion_default = Quaternion( (0,0,0,0) )[:]

class SvListInputQuaternionEntry(bpy.types.PropertyGroup):

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value

    def get_elem_SCALARVECTOR(self):
        #q = Quaternion(self.elem)
        return self.elem

    def set_elem_SCALARVECTOR(self, value):
        self.elem = value

    def get_elem_EULER(self):
        q = Quaternion(self.elem)
        angle = q.to_euler('XYZ')
        #return self.elem
        return angle

    def set_elem_EULER(self, value):
        euler = Euler( value, 'XYZ' )
        q = euler.to_quaternion()
        self.elem = q[:]
        pass

    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    quaternion_euler_order: EnumProperty(
        name="Euler Order",
        description="Order of the Euler rotations",
        default="XYZ",
        items=quaternion_euler_orders,
    ) # type: ignore
    
    elem: FloatVectorProperty(
        name = "Vector",
        default = quaternion_default, #(0., 0., 0., 0., ),
        size=4,
        update=update_entry,
    ) # type: ignore

    NONE: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0., 0., ),
        subtype='NONE',
        size=4,
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore

    
    WXYZ: FloatVectorProperty(
        name = "Float number",
        default = (0., 0., 0., 0., ),
        subtype='QUATERNION',
        size=4,
        #update=update_entry,
    ) # type: ignore

    def get_WXYZ(self):
        return self.WXYZ

    def set_WXYZ(self, value):
        self.WXYZ = value
        q = Quaternion(value)
        euler = q.to_euler('XYZ')
        self.elem = q[:]
        self.SCALARVECTOR = q[:]
        self.EULER = euler
        self.AXISANGLE = q[:]
        pass

    WXYZ_UI: FloatVectorProperty(
        name = "Quaternion",
        default = (0., 0., 0., 0., ),
        subtype='QUATERNION',
        size=4,
        get = get_WXYZ, set = set_WXYZ,
        update=update_entry,
    ) # type: ignore

    EULER: FloatVectorProperty(
        name = "Eular Angle (XYZ)",
        default = (0., 0., 0., ),
        subtype='EULER',
        size=3,
    ) # type: ignore

    def get_EULER(self):
        return self.EULER

    def set_EULER(self, value):
        self.EULER = value
        euler = Euler( value, 'XYZ' )
        q = euler.to_quaternion()
        self.elem = q[:]
        self.WXYZ = q[:]
        self.SCALARVECTOR = q[:]
        self.AXISANGLE = q[:]
        pass

    EULER_UI: FloatVectorProperty(
        name = "Eular Angle (XYZ)",
        default = (0., 0., 0., ),
        subtype='EULER',
        size=3,
        get = get_EULER, set = set_EULER,
        update=update_entry,
    ) # type: ignore

    SCALARVECTOR: FloatVectorProperty(
        name = "Scalar Vector (Scalar, XYZ)",
        default = (0., 0., 0., 0, ),
        subtype='QUATERNION',
        size=4,
    ) # type: ignore

    def get_SCALARVECTOR(self):
        return self.SCALARVECTOR

    def set_SCALARVECTOR(self, value):
        self.SCALARVECTOR = value
        q = Quaternion(value)
        euler = q.to_euler('XYZ')
        self.elem = q[:]
        self.WXYZ = q[:]
        self.EULER = euler
        self.AXISANGLE = q[:]
        pass

    SCALARVECTOR_UI: FloatVectorProperty(
        name = "Scalar Vector",
        default = (0., 0., 0., 0., ),
        subtype='QUATERNION',
        size=4,
        get = get_SCALARVECTOR, set = set_SCALARVECTOR,
        update=update_entry,
    ) # type: ignore
    

    AXISANGLE: FloatVectorProperty(
        name = "Angle-Axis (Angle, XYZ)",
        default = (0., 0., 0., 0, ),
        subtype='QUATERNION',
        size=4,
    ) # type: ignore

    def get_AXISANGLE(self):
        return self.AXISANGLE

    def set_AXISANGLE(self, value):
        self.AXISANGLE = value
        q = Quaternion( value[1:], value[1] )
        euler = q.to_euler('XYZ')
        self.elem = q[:]
        self.WXYZ = q[:]
        self.SCALARVECTOR = q[:]
        self.EULER = euler
        pass

    AXISANGLE_UI: FloatVectorProperty(
        name = "Angle-Axis (Angle, XYZ)",
        default = (0., 0., 0., 0., ),
        subtype='AXISANGLE',
        size=4,
        get = get_AXISANGLE, set = set_AXISANGLE,
        update=update_entry,
    ) # type: ignore

class SvListInputColorEntry(bpy.types.PropertyGroup):

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value

    def get_elem_without_alpha(self):
        return self.elem[:3]

    def set_elem_without_alpha(self, value):
        self.elem = (*value[:3], self.elem[3] )

    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: FloatVectorProperty(
        name = "Vector",
        default = (1., 1., 1., 1.),
        size=4, min=0., max=1.,
        update=update_entry,
    ) # type: ignore

    elem_WITHOUT_ALPHA: FloatVectorProperty(
        name = "Vector",
        default = (1., 1., 1., ),
        size=3, min=0., max=1.,
        update=update_entry,
    ) # type: ignore

    NONE: FloatVectorProperty(
        name = "Color as vector",
        default = (1., 1., 1., 1.),
        subtype='NONE',
        size=4, min=0., max=1.,
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    NONE_WITHOUT_ALPHA: FloatVectorProperty(
        name = "Color as vector",
        default = (1., 1., 1., ),
        subtype='NONE',
        size=3, min=0., max=1.,
        get = get_elem_without_alpha, set = set_elem_without_alpha,
        update=update_entry,
    ) # type: ignore

    COLOR: FloatVectorProperty(
        name = "Color",
        default = (1., 1., 1., 1.),
        size=4, min=0., max=1.,
        subtype='COLOR',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    COLOR_WITHOUT_ALPHA: FloatVectorProperty(
        name = "Color",
        default = (1., 1., 1., ),
        size=3, min=0., max=1.,
        subtype='COLOR',
        get = get_elem_without_alpha, set = set_elem_without_alpha,
        update=update_entry,
    ) # type: ignore

    COLOR_GAMMA: FloatVectorProperty(
        name = "Color",
        default = (1., 1., 1., 1.),
        size=4, min=0., max=1.,
        subtype='COLOR_GAMMA',
        get = get_elem, set = set_elem,
        update=update_entry,
    ) # type: ignore

    COLOR_GAMMA_WITHOUT_ALPHA: FloatVectorProperty(
        name = "Color Gamma without alpha",
        default = (1., 1., 1., ),
        size=3, min=0., max=1.,
        subtype='COLOR_GAMMA',
        get = get_elem_without_alpha, set = set_elem_without_alpha,
        update=update_entry,
    ) # type: ignore


class SvListInputStringEntry(bpy.types.PropertyGroup):

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value

    def update_entry(self, context):
        # TODO: Если нод не активен, то изменения текста через контролы выбора файла или каталога
        # не работают и выходные значения нода не обновляются. Как-то надо найти нод, к
        # которому относятся эти контролы.
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        elif hasattr(context, 'active_node'):
            updateNode(context.active_node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
        
    elem: StringProperty(
        name = "String of text",
        default = "",
        update=update_entry,
    ) # type: ignore
    
    NONE: StringProperty(
        name = "String of text",
        default = "",
        subtype='NONE',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore
    
    FILE_PATH: StringProperty(
        name = "String of text",
        default = "",
        subtype='FILE_PATH',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore
    
    DIR_PATH: StringProperty(
        name = "String of text",
        default = "",
        subtype='DIR_PATH',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore
    
    FILE_NAME: StringProperty(
        name = "String of text",
        default = "",
        subtype='FILE_NAME',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore
    
    BYTE_STRING: StringProperty(
        name = "String of text",
        default = "",
        subtype='BYTE_STRING',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore
    
    PASSWORD: StringProperty(
        name = "String of text",
        default = "",
        subtype='PASSWORD',
        get = get_elem, set = set_elem,
        #update=update_entry,
    ) # type: ignore

class SvListInputMatrixEntry(bpy.types.PropertyGroup):

    id_matrix = (1.0, 0.0, 0.0, 0.0,
                 0.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 0.0,
                 0.0, 0.0, 0.0, 1.0)

    # share data with elements
    def get_elem(self):
        return self.elem

    def set_elem(self, value):
        self.elem = value

    def update_entry(self, context):
        if hasattr(context, 'node'):
            updateNode(context.node, context)
        else:
            sv_logger.debug("Node is not defined in this context, so will not update the node.")

    item_enable : BoolProperty(
        name = "Mask",
        description = "On - Add Element in the output socket,\nOff - Do not add Element in the output socket",
        default=True,
        update=update_entry,
        ) # type: ignore
    
    elem: FloatVectorProperty(
        name = "Vector",
        default = id_matrix,
        #subtype='MATRIX',
        size=16,
        update=update_entry,
    ) # type: ignore

    MATRIX: FloatVectorProperty(
        name = "Matrix elem",
        default = id_matrix,
        size=16,
        #subtype='MATRIX',
        update=update_entry,
    ) # type: ignore

    # share data with elements
    def get_MATRIX(self):
        return self.MATRIX

    def set_MATRIX(self, value):
        mat = np.array(value).reshape(-1,4).tolist()
        # for i in range(4):
        #     arr_i = []
        #     mat.append(arr_i)
        #     for j in range(i, 16, 4):
        #         arr_i.append(value[j])

        self.MATRIX = value
        mat = Matrix(mat)
        T, R, S = mat.decompose()
        
        self.EULER_LOCATION[:] = T[:]
        self.EULER_ANGLE[:]    = R.to_euler("XYZ")[:]
        self.EULER_SCALE[:]    = S[:]

        self.AXISANGLE_LOCATION[:] = T[:]
        self.AXISANGLE_SCALE[:]    = S[:]
        self.AXISANGLE_VECTOR[:]   = R.axis[:]
        self.AXISANGLE_ANGLE       = R.angle

        pass

    MATRIX_UI: FloatVectorProperty(
        name = "Matrix elem",
        default = id_matrix,
        size=16,
        precision=3,
        #subtype='MATRIX',
        get = get_MATRIX, set = set_MATRIX,
        update=update_entry,
    ) # type: ignore

    ####### EULER_LOCATION, euler_scale, euler_angle_xyz #####################################################################################
    EULER_LOCATION: FloatVectorProperty(
        name = "Eular Angle (XYZ)",
        default = (0., 0., 0., ),
        #subtype='EULER',
        size=3,
    ) # type: ignore

    def update_by_EULER(self):
        # translation
        mat_t = Matrix().Identity(4)
        mat_t[0][3] = self.EULER_LOCATION[0]
        mat_t[1][3] = self.EULER_LOCATION[1]
        mat_t[2][3] = self.EULER_LOCATION[2]

        # rotation
        mat_r = Euler(self.EULER_ANGLE[:], "XYZ").to_matrix().to_4x4()
        # scale
        mat_s = Matrix().Identity(4)
        mat_s[0][0] = self.EULER_SCALE[0]
        mat_s[1][1] = self.EULER_SCALE[1]
        mat_s[2][2] = self.EULER_SCALE[2]
        # composite matrix
        mat = mat_t @ mat_r @ mat_s
        T, R, S = mat.decompose()

        self.MATRIX = list(chain( *mat.row ))
        # skip EULER
        self.AXISANGLE_LOCATION[:] = T[:]
        self.AXISANGLE_SCALE[:]    = S[:]
        self.AXISANGLE_VECTOR[:]   = R.axis[:]
        self.AXISANGLE_ANGLE       = R.angle

    def get_EULER_LOCATION(self):
        return self.EULER_LOCATION

    def set_EULER_LOCATION(self, value):
        self.EULER_LOCATION = value
        self.update_by_EULER()
        pass

    EULER_LOCATION_UI_NONE: FloatVectorProperty(
        name = "Eular Location (XYZ)",
        default = (0., 0., 0., ),
        precision=3,
        #subtype='NONE',
        size=3,
        get = get_EULER_LOCATION, set = set_EULER_LOCATION,
        update=update_entry,
    ) # type: ignore
    EULER_LOCATION_UI_TRANSLATION: FloatVectorProperty(
        name = "Eular Location (XYZ)",
        default = (0., 0., 0., ),
        subtype='TRANSLATION',
        precision=3,
        size=3,
        get = get_EULER_LOCATION, set = set_EULER_LOCATION,
        update=update_entry,
    ) # type: ignore

    ####### /EULER_LOCATION, euler_scale, euler_angle_xyz #####################################################################################
    ####### euler_location, EULER_SCALE, euler_angle_xyz #####################################################################################

    EULER_SCALE: FloatVectorProperty(
        name = "Eular Scale (XYZ)",
        default = (1., 1., 1., ),
        #subtype='EULER',
        size=3,
    ) # type: ignore

    def get_EULER_SCALE(self):
        return self.EULER_SCALE

    def set_EULER_SCALE(self, value):
        self.EULER_SCALE = value
        self.update_by_EULER()
        pass

    EULER_SCALE_UI: FloatVectorProperty(
        name = "Eular Scale (XYZ)",
        default = (1., 1., 1., ),
        precision=3,
        #subtype='EULER',
        size=3,
        get = get_EULER_SCALE, set = set_EULER_SCALE,
        update=update_entry,
    ) # type: ignore
    
    ####### euler_location, /EULER_SCALE, euler_angle_xyz #####################################################################################
    ####### euler_location,  euler_scale, EULER_ANGLE_XYZ #####################################################################################

    EULER_ANGLE: FloatVectorProperty(
        name = "Eular Angle (XYZ)",
        default = (0., 0., 0., ),
        subtype='EULER',
        size=3,
    ) # type: ignore

    def get_EULER_ANGLE(self):
        return self.EULER_ANGLE

    def set_EULER_ANGLE(self, value):
        self.EULER_ANGLE = value
        self.update_by_EULER()
        pass

    EULER_ANGLE_UI: FloatVectorProperty(
        name = "Eular Scale (XYZ)",
        default = (0., 0., 0., ),
        subtype='EULER',
        size=3,
        get = get_EULER_ANGLE, set = set_EULER_ANGLE,
        update=update_entry,
    ) # type: ignore
    ####### euler_location,  euler_scale, /EULER_ANGLE_XYZ #####################################################################################

    ####### AXISANGLE_LOCATION, axisangle_scale, axisangle_vector, axisangle_angle #########################################################################
    AXISANGLE_LOCATION: FloatVectorProperty(
        name = "Axis Angle Location (XYZ)",
        default = (0., 0., 0., ),
        #subtype='EULER',
        size=3,
    ) # type: ignore

    def update_by_AXISANGLE(self):
        # translation
        mat_t = Matrix().Identity(4)
        mat_t[0][3] = self.AXISANGLE_LOCATION[0]
        mat_t[1][3] = self.AXISANGLE_LOCATION[1]
        mat_t[2][3] = self.AXISANGLE_LOCATION[2]

        # rotation
        mat_r = Quaternion(self.AXISANGLE_VECTOR[:], self.AXISANGLE_ANGLE).to_matrix().to_4x4()
        # scale
        mat_s = Matrix().Identity(4)
        mat_s[0][0] = self.AXISANGLE_SCALE[0]
        mat_s[1][1] = self.AXISANGLE_SCALE[1]
        mat_s[2][2] = self.AXISANGLE_SCALE[2]
        # composite matrix
        mat = mat_t @ mat_r @ mat_s
        T, R, S = mat.decompose()

        self.MATRIX = list(chain( *mat.row ))
        self.EULER_LOCATION[:] = T[:]
        self.EULER_ANGLE[:]    = R.to_euler("XYZ")[:]
        self.EULER_SCALE[:]    = S[:]
        # skip AXISANGLE
        pass

    def get_AXISANGLE_LOCATION(self):
        return self.AXISANGLE_LOCATION

    def set_AXISANGLE_LOCATION(self, value):
        self.AXISANGLE_LOCATION = value
        self.update_by_AXISANGLE()
        pass

    AXISANGLE_LOCATION_UI_NONE: FloatVectorProperty(
        name = "Axis-Angle Location (XYZ)",
        default = (0., 0., 0., ),
        #subtype='EULER',
        size=3,
        get = get_AXISANGLE_LOCATION, set = set_AXISANGLE_LOCATION,
        update=update_entry,
    ) # type: ignore
    AXISANGLE_LOCATION_UI_TRANSLATION: FloatVectorProperty(
        name = "Axis-Angle Location (XYZ)",
        default = (0., 0., 0., ),
        subtype='TRANSLATION',
        size=3,
        get = get_AXISANGLE_LOCATION, set = set_AXISANGLE_LOCATION,
        update=update_entry,
    ) # type: ignore

    ####### /AXISANGLE_LOCATION, axisangle_scale, axisangle_vector, axisangle_angle #########################################################################
    #######  axisangle_location, AXISANGLE_SCALE, axisangle_vector, axisangle_angle #########################################################################
    AXISANGLE_SCALE: FloatVectorProperty(
        name = "Axis-Angle Scale (XYZ)",
        default = (1., 1., 1., ),
        #subtype='EULER',
        size=3,
    ) # type: ignore

    def get_AXISANGLE_SCALE(self):
        return self.AXISANGLE_SCALE

    def set_AXISANGLE_SCALE(self, value):
        self.AXISANGLE_SCALE = value
        self.update_by_AXISANGLE()
        pass

    AXISANGLE_SCALE_UI: FloatVectorProperty(
        name = "Eular Scale (XYZ)",
        default = (1., 1., 1., ),
        #subtype='EULER',
        size=3,
        get = get_AXISANGLE_SCALE, set = set_AXISANGLE_SCALE,
        update=update_entry,
    ) # type: ignore

    ####### axisangle_location, /AXISANGLE_SCALE, axisangle_vector, axisangle_angle #########################################################################
    ####### axisangle_location,  axisangle_scale, AXISANGLE_VECTOR, axisangle_angle #########################################################################
    AXISANGLE_VECTOR: FloatVectorProperty(
        name = "Axis (XYZ)",
        default = (0., 0., 1., ),
        #subtype='EULER',
        size=3,
    ) # type: ignore

    def get_AXISANGLE_VECTOR(self):
        return self.AXISANGLE_VECTOR

    def set_AXISANGLE_VECTOR(self, value):
        self.AXISANGLE_VECTOR = value
        self.update_by_AXISANGLE()
        pass

    AXISANGLE_VECTOR_UI: FloatVectorProperty(
        name = "Eular Angle (XYZ)",
        default = (0., 0., 1., ),
        #subtype='EULER',
        size=3,
        get = get_AXISANGLE_VECTOR, set = set_AXISANGLE_VECTOR,
        update=update_entry,
    ) # type: ignore

    ####### axisangle_location, axisangle_scale, /AXISANGLE_VECTOR, axisangle_angle #########################################################################
    ####### axisangle_location, axisangle_scale,  axisangle_vector, AXISANGLE_ANGLE #########################################################################
    AXISANGLE_ANGLE: FloatProperty(
        name = "Axis-Angle Angle",
        default = 0.,
    ) # type: ignore

    def get_AXISANGLE_ANGLE(self):
        return self.AXISANGLE_ANGLE

    def set_AXISANGLE_ANGLE(self, value):
        self.AXISANGLE_ANGLE = value
        self.update_by_AXISANGLE()
        pass

    AXISANGLE_ANGLE_UI: FloatProperty(
        name = "Axis-Angle Angle",
        default = 0.,
        subtype='ANGLE',
        get = get_AXISANGLE_ANGLE, set = set_AXISANGLE_ANGLE,
        update=update_entry,
    ) # type: ignore

    ####### axis_angle_location, /AXISANGLE_SCALE, axisangle_vector, AXISANGLE_ANGLE #########################################################################


class SvCopyTextToClipboard(bpy.types.Operator, SvGenericNodeLocator):
    ''' Copy node's data by data mode (bool, int, float, vector, quaternion, color, strings) to the clipboard '''
    bl_idname = "node.sverchok_copy_text_to_clipboard"
    bl_label = "sverchok: copy data to the clipboard as text"
    # bl_options = {'REGISTER', 'UNDO'}

    # def execute(self, context):
    #     node = self.get_node(context)
    #     context.window_manager.clipboard = self.text
    #     return {'FINISHED'}
    
    def sv_execute(self, context, node):
        if hasattr(node, 'dataAsString')==True:
            text = node.dataAsString()
            context.window_manager.clipboard = text
            ShowMessageBox("Data copied to the Clipboard")
        pass

class SvPasteTextFromClipboard(bpy.types.Operator, SvGenericNodeLocator):
    '''Paste data into the node by mode (bool, int, float, vector, quaternion, color, strings) from the clipboard text'''
    bl_idname = "node.sverchok_paste_text_from_clipboard"
    bl_label = "sverchok: paste data from the clipboard"
    bl_options = {'REGISTER', 'UNDO'}

    # def execute(self, context):
    #     node = self.get_node(context)
    #     context.window_manager.clipboard = self.text
    #     return {'FINISHED'}
    
    def sv_execute(self, context, node):
        if hasattr(node, 'convertTextToData')==True:
            text = context.window_manager.clipboard
            node.convertTextToData(text)
        pass

class SvUpdateTextInListInputNode(bpy.types.Operator, SvGenericNodeLocator):
    '''Update text data. This is a glitch on Blender. If node is not active then edited text fields are not updated automatically'''
    bl_idname = "node.sverchok_update_text_in_list_input_node"
    bl_label = "sverchok: update text in list input node"
    # bl_options = {'REGISTER', 'UNDO'}

    def sv_execute(self, context, node):
        if hasattr(node, 'updateTextData')==True:
            node.updateTextData(context)
        pass


def Correct_ListInput_Length(self, context):
    if(self.mode=='BOOL_LIST_MODE'):
        correct_collection_length(self.bool_list_items, self.bool_list_counter)
    elif(self.mode=='INT_LIST_MODE'):
        correct_collection_length(self.int_list_items, self.int_list_counter)
    elif(self.mode=='FLOAT_LIST_MODE'):
        correct_collection_length(self.float_list_items, self.float_list_counter)
    elif(self.mode=='VECTOR_LIST_MODE'):
        correct_collection_length(self.vector_list_items, self.vector_list_counter)
    elif(self.mode=='QUATERNION_LIST_MODE'):
        correct_collection_length(self.quaternion_list_items, self.quaternion_list_counter)
    elif(self.mode=='MATRIX_LIST_MODE'):
        correct_collection_length(self.matrix_list_items, self.matrix_list_counter)
    elif(self.mode=='COLOR_LIST_MODE'):
        correct_collection_length(self.color_list_items, self.color_list_counter)
    elif(self.mode=='STRING_LIST_MODE'):
        correct_collection_length(self.string_list_items, self.string_list_counter)
    else:
        raise Exception(f"[func: Correct_ListInput_Length] unknown mode {self.mode}.")
    updateNode(self, context)
    pass

def invert_ListInput_Mask(self, context):
    if(self.mode=='BOOL_LIST_MODE'):
        items = self.bool_list_items
    elif(self.mode=='INT_LIST_MODE'):
        items = self.int_list_items
    elif(self.mode=='FLOAT_LIST_MODE'):
        items = self.float_list_items
    elif(self.mode=='VECTOR_LIST_MODE'):
        items = self.vector_list_items
    elif(self.mode=='QUATERNION_LIST_MODE'):
        items = self.quaternion_list_items
    elif(self.mode=='MATRIX_LIST_MODE'):
        items = self.matrix_list_items
    elif(self.mode=='COLOR_LIST_MODE'):
        items = self.color_list_items
    elif(self.mode=='STRING_LIST_MODE'):
        items = self.string_list_items
    else:
        raise Exception(f"[func: invert_ListInput_Mask] unknown mode {self.mode}.")
    for item in items:
        item.item_enable = not(item.item_enable)
    updateNode(self, context)
    pass

# def copy_to_clipboard(self, context):
#     #context.window_manager.clipboard = '1234'
#     if self.mode == 'BOOL_LIST_MODE':
#         lst = [elem.elem for elem in self.bool_list_items if elem.item_enable==True]
#     elif self.mode == 'INT_LIST_MODE':
#         lst = [elem.elem for elem in self.int_list_items if elem.item_enable==True]
#     elif self.mode == 'FLOAT_LIST_MODE':
#         lst = [elem.elem for elem in self.float_list_items if elem.item_enable==True]
#     elif self.mode == 'VECTOR_LIST_MODE':
#         lst = []
#         [lst.extend(tuple(Vector(elem.elem))) for elem in self.color_list_items if elem.item_enable==True]
#         #lst = [str(tuple(Vector(elem.elem))) for elem in self.vector_list_items if elem.item_enable==True]
#     elif self.mode == 'COLOR_LIST_MODE':
#         lst = []
#         [lst.extend(tuple(elem.elem)) for elem in self.color_list_items if elem.item_enable==True]
#     else:
#         raise Exception(f"[func: process] unknown mode {self.mode}.")
#     #text = ','.join([str(elem) for elem in lst])
#     text = ', '.join([str(e) for e in lst])
#     context.window_manager.clipboard = text
#     pass

class SvListInputNodeMK2(Show3DProperties, SverchCustomTreeNode, bpy.types.Node):
    ''' Create a flat list of Integers, Floats, or Vectors.
    int: [[1,5,42]],
    float: [[1.0, -5.42, 12.0]]
    vectors: [[ (0,0,0), (1,2,-2)]]
    '''
    bl_idname = 'SvListInputNodeMK2'
    bl_label = 'List Input'
    bl_icon = 'OUTLINER_OB_EMPTY'
    sv_icon = 'SV_LIST_INPUT'

    def wrapper_tracked_ui_draw_op(self, layout_element, operator_idname, **keywords):
        """
        this wrapper allows you to track the origin of a clicked operator, by automatically passing
        the node_name and tree_name to the operator.

        example usage:

            row.separator()
            self.wrapper_tracked_ui_draw_op(row, "node.view3d_align_from", icon='CURSOR', text='')

        """
        op = layout_element.operator(operator_idname, **keywords)
        op.node_name = self.name
        op.tree_name = self.id_data.name
        return op


    bool_list_counter: IntProperty(
        name='bool_list_counter',
        description='boolean',
        default=True,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    bool_list_items : CollectionProperty(type=SvListInputBoolEntry) # type: ignore

    int_list_counter: IntProperty(
        name='int_list_counter',
        description='integer number',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    int_list_items : CollectionProperty(type=SvListInputIntEntry) # type: ignore

    float_list_counter: IntProperty(
        name='float_list_counter',
        description='integer number',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    float_list_items : CollectionProperty(type=SvListInputFloatEntry) # type: ignore

    vector_list_counter: IntProperty(
        name='vector_list_counter',
        description='vectors',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    vector_list_items : CollectionProperty(type=SvListInputVectorEntry) # type: ignore

    matrix_list_counter: IntProperty(
        name='matrix_list_counter',
        description='matrixes',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    matrix_list_items : CollectionProperty(type=SvListInputMatrixEntry) # type: ignore

    quaternion_list_counter: IntProperty(
        name='quaternion_list_counter',
        description='vectors',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    quaternion_list_items : CollectionProperty(type=SvListInputQuaternionEntry) # type: ignore

    color_list_counter: IntProperty(
        name='color_list_counter',
        description='Colors',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    color_list_items : CollectionProperty(type=SvListInputColorEntry) # type: ignore

    string_list_counter: IntProperty(
        name='string_list_counter',
        description='String',
        default=1,
        min=1,
        update=Correct_ListInput_Length)  # type: ignore
    string_list_items : CollectionProperty(type=SvListInputStringEntry) # type: ignore

    invert_mask: BoolProperty(
        name = "Mask",
        description = "Invert list of mask. Disabled if input socket 'Mask' are connected",
        default=True,
        update=invert_ListInput_Mask,
    ) # type: ignore

    # copy_clipboard: BoolProperty(
    #     name = "Mask",
    #     description = "Invert list of mask",
    #     default=True,
    #     update=copy_to_clipboard,
    # ) # type: ignore

    
    mask_modes = [
            ('BOOLEAN', "Booleans", "Boolean values (0/1) as mask [[0,1,0,0,1,1],[1,1,0,0,1],...]. Has no influence if socket is not connected (All data are used)", 0),
            ('INDEXES', "Indexes", "Indexes as mask [[1,2,0,4],[0,1,4,5,7],..]. Has no influence if socket is not connected (All data are used)", 1),
        ]
    mask_mode : EnumProperty(
        name = "Mask mode",
        items = mask_modes,
        default = 'BOOLEAN',
        #update = updateMaskMode
        update = updateNode
        ) # type: ignore

    mask_inversion : BoolProperty(
        name = "Invert",
        default = False,
        description="Invert mask. Has no influence if socket is not connected (All data are used)",
        update = updateNode) # type: ignore
    
    def draw_mask_in_socket(self, socket, context, layout):
        grid = layout.grid_flow(row_major=True, columns=2)
        if not socket.is_linked:
            grid.enabled = False
        col2 = grid.column()
        col2_row1 = col2.row()
        col2_row1.alignment='LEFT'
        if socket.is_linked:
            col2_row1.label(text=f"Mask. {socket.objects_number or ''}:")
        else:
            col2_row1.label(text=f"Mask: (not connected)")
        col2_row2 = col2.row()
        col2_row2.alignment='LEFT'
        col2_row2.column(align=True).prop(self, "mask_inversion")
        col3 = grid.column()
        col3.prop(self, "mask_mode", expand=True)

    def changeMode(self, context):
        if self.mode == 'BOOL_LIST_MODE':
            self.set_output_socketype(['SvStringsSocket'])
            self.outputs['data_output'].label = 'Boolean'
        elif self.mode == 'INT_LIST_MODE':
            self.set_output_socketype(['SvStringsSocket'])
            self.outputs['data_output'].label = 'Integers'
        elif self.mode == 'FLOAT_LIST_MODE':
            self.set_output_socketype(['SvStringsSocket'])
            self.outputs['data_output'].label = 'Floats'
        elif self.mode == 'VECTOR_LIST_MODE':
            self.set_output_socketype(['SvVerticesSocket'])
            self.outputs['data_output'].label = 'Vectors'
        elif self.mode == 'QUATERNION_LIST_MODE':
            self.set_output_socketype(['SvQuaternionSocket'])
            self.outputs['data_output'].label = 'Quaternions'
        elif self.mode == 'MATRIX_LIST_MODE':
            self.set_output_socketype(['SvMatrixSocket'])
            self.outputs['data_output'].label = 'Matrixes'
        elif self.mode == 'COLOR_LIST_MODE':
            self.set_output_socketype(['SvColorSocket'])
            self.outputs['data_output'].label = 'Colors'
        elif self.mode == 'STRING_LIST_MODE':
            self.set_output_socketype(['SvTextSocket'])
            self.outputs['data_output'].label = 'Strings'
        else:
            raise Exception(f"[func: changeMode] unknown mode {self.mode}.")
        Correct_ListInput_Length(self, context)
        pass

    modes = [
        ("BOOL_LIST_MODE", "Bool", "Boolean", "IMAGE_ALPHA", 0),
        ("INT_LIST_MODE", "Int", "Integer", "IPO_CONSTANT", 1),
        ("FLOAT_LIST_MODE", "Float", "Float", "IPO_LINEAR", 2),
        ("VECTOR_LIST_MODE", "Vector", "Vector", "ORIENTATION_GLOBAL", 3),
        ("QUATERNION_LIST_MODE", "Quaternion", "Quaternion", "CURVE_PATH", 4),
        ("MATRIX_LIST_MODE", "Matrix", "Matrix", "MOD_LATTICE", 5),
        ("COLOR_LIST_MODE", "Color", "Color", "COLOR", 6),
        ("STRING_LIST_MODE", "Text", "Text", "SORTALPHA", 7),
    ]

    mode: EnumProperty(
        description="Data types of node",
        items=modes,
        default='INT_LIST_MODE',
        update=changeMode
    ) # type: ignore

    subtypes_vector = [
        # # https://docs.blender.org/api/current/bpy_types_enum_items/property_subtype_number_array_items.html#rna-enum-property-subtype-number-array-items
        ("NONE", "None", "None", "", 0),
        #("COLOR", "Color", "Color", "", 1),
        ("TRANSLATION", "Translation", "Translation", "", 2),
        ("DIRECTION", "Direction", "Direction", "", 3),
        ("VELOCITY", "Velocity", "Velocity", "", 4),
        ("ACCELERATION", "Acceleration", "Acceleration", "", 5),
        #("MATRIX", "Matrix", "Matrix", "", 6),
        ("EULER", "Euler Angles XYZ", "Euler Angles", "", 7),
        #("QUATERNION", "Quaternion", "Quaternion", "", 8), - moved to mode
        ("AXISANGLE", "Axis-Angle", "Axis-Angle", "", 9),
        ("XYZ", "XYZ", "XYZ", "", 10),
        ("XYZ_LENGTH", "XYZ Length", "XYZ Length", "", 11),
        #("COLOR_GAMMA", "Color Gamma", "Color Gamma", "", 12),
        ("COORDINATES", "Coordinates", "Coordinates", "", 13),
        # ("LAYER", "Layer", "Layer", "", 14),
        # ("LAYER_MEMBER", "Layer Member", "Layer Member", "", 15),
    ]

    subtype_vector: EnumProperty(
        name='Subtypes',
        items=subtypes_vector,
        default='NONE',
        #update=changeSubType
    ) # type: ignore

    subtypes_color = [
        # # https://docs.blender.org/api/current/bpy_types_enum_items/property_subtype_number_array_items.html#rna-enum-property-subtype-number-array-items
        ("NONE", "None", "None (as vector)", "CHECKBOX_DEHLT", 0),
        ("COLOR", "Color", "Color", "RESTRICT_COLOR_ON", 1),
        ("COLOR_GAMMA", "Color Gamma", "Color Gamma", "COLOR", 12),
    ]

    subtype_color: EnumProperty(
        name='Subtypes',
        items=subtypes_color,
        default='COLOR',
        #update=changeSubType
    ) # type: ignore

    subtypes_float = [
        # https://docs.blender.org/api/current/bpy_types_enum_items/property_subtype_number_items.html#rna-enum-property-subtype-number-items
        ("NONE", "None", "None", "", 0),
        #("PERCENTAGE", "Percentage", "Percentage", "", 1),
        #("FACTOR", "Factor", "Factor", "", 2),
        ("ANGLE", "Angle", "Angle", "", 3),
        #("TIME", "Time", "Time", "", 4),
        ("TIME_ABSOLUTE", "Time(abs)", "Time absolute", "", 5),
        ("DISTANCE", "Distance", "Distance", "", 6),
        ("DISTANCE_CAMERA", "Distance(cam)", "Distance camera", "", 7),
        ("TEMPERATURE", "Temperature", "Temperature", "", 8),
        #("WAVELENGTH", "Wavelength", "Wavelength", "", 9), - Error in Blender ?
        ("POWER", "Power", "Power", "", 10),
    ]

    subtype_float: EnumProperty(
        name='Subtypes',
        items=subtypes_float,
        default='NONE',
        #update=changeSubTypeFloat
    ) # type:  ignore


    subtypes_string = [
        # https://docs.blender.org/api/current/bpy_types_enum_items/property_subtype_string_items.html#rna-enum-property-subtype-string-items
        ("NONE", "None", "None", "", 0),
        ("FILE_PATH", "File Path", "File Path", "", 1),
        ("DIR_PATH", "Directory Path", "Directory Path", "", 2),
        ("FILE_NAME", "File Name", "File Name", "", 3),
        ("BYTE_STRING", "Byte String", "Byte String", "", 4),
        ("PASSWORD", "Password", "Password", "A string that is displayed hidden (’********’)", 5),
    ]
    
    subtype_string: EnumProperty(
        name='Subtypes',
        items=subtypes_string,
        default='NONE',
        #update=changeSubTypeFloat
    ) # type:  ignore

    # share data with elements
    def get_unit_system(self):
        system_unit = bpy.data.scenes['Scene'].unit_settings.system
        res = [e for e in self.unit_systems if e[0]==system_unit][0][4]
        return res

    def set_unit_system(self, value):
        bpy.data.scenes['Scene'].unit_settings.system = self.unit_systems[value][0]
        pass

    unit_systems = [
        ("NONE", "None", "None", "", 0),
        ("METRIC", "Metric", "Metric", "", 1),
        ("IMPERIAL", "Imperial", "Imperial]", "", 2),
    ]
    unit_system: EnumProperty(
        name='Unit System',
        items=unit_systems,
        default='NONE',
        get=get_unit_system, set=set_unit_system,
        #update=changeSubTypeFloat
    ) # type:  ignore

    length_units_metric = [
        ("ADAPTIVE", "Adaptive", "Adaptive", "", 0),
        ("KILOMETERS", "Kilometers", "Kilometers", "10^3", 1),
        ("METERS", "Meters", "Meters", "1", 2),
        ("CENTIMETERS", "Centimeters", "Centimeters", "10^-2", 3),
        ("MILLIMETERS", "Millimeters", "Millimeters", "", 4),
        ("MICROMETERS", "Micrometers", "Micrometers", "", 5),
    ]

    def get_length_unit_metric(self):
        length_unit = bpy.data.scenes['Scene'].unit_settings.length_unit
        res = [e for e in self.length_units_metric if e[0]==length_unit][0][4]
        return res

    def set_length_unit_metric(self, value):
        bpy.data.scenes['Scene'].unit_settings.length_unit = self.length_units_metric[value][0]
        pass

    length_unit_metric: EnumProperty(
        name='Length unit metric',
        items=length_units_metric,
        default='METERS',
        get=get_length_unit_metric, set=set_length_unit_metric,
        #update=changeSubTypeFloat
    ) # type:  ignore

    def get_length_unit_imperial(self):
        length_unit = bpy.data.scenes['Scene'].unit_settings.length_unit
        res = [e for e in self.length_units_imperial if e[0]==length_unit][0][4]
        return res

    def set_length_unit_imperial(self, value):
        bpy.data.scenes['Scene'].unit_settings.length_unit = self.length_units_imperial[value][0]
        pass

    length_units_imperial = [
        ("ADAPTIVE", "Adaptive", "Adaptive", "", 0),
        ("MILES", "Miles", "Miles", "", 1),
        ("FEET", "Feet", "Feet", "1", 2),
        ("INCHES", "Inches", "Inches", "", 3),
        ("THOU", "Thou", "Thow", "", 4),
    ]
    length_unit_imperial: EnumProperty(
        name='Length unit imperial',
        items=length_units_imperial,
        default='FEET',
        get=get_length_unit_imperial, set=set_length_unit_imperial,
        #update=changeSubTypeFloat
    ) # type:  ignore

    def get_rotation_unit(self):
        system_rotation = bpy.data.scenes['Scene'].unit_settings.system_rotation
        res = [e for e in self.system_rotations if e[0]==system_rotation][0][4]
        return res

    def set_rotation_unit(self, value):
        bpy.data.scenes['Scene'].unit_settings.system_rotation = self.system_rotations[value][0]
        pass

    system_rotations = [
        ("DEGREES", "Degrees", "Degrees", "", 0),
        ("RADIANS", "Radians", "Radians", "", 1),
    ]
    system_rotation: EnumProperty(
        name='Rotation unit',
        items=system_rotations,
        default='RADIANS',
        get=get_rotation_unit, set=set_rotation_unit,
        #update=changeSubTypeFloat
    ) # type:  ignore

    def get_temperature_unit_metric(self):
        temperature_unit = bpy.data.scenes['Scene'].unit_settings.temperature_unit
        res = [e for e in self.temperature_units_metric if e[0]==temperature_unit][0][4]
        return res

    def set_temperature_unit_metric(self, value):
        bpy.data.scenes['Scene'].unit_settings.temperature_unit = self.temperature_units_metric[value][0]
        pass

    temperature_units_metric = [
        ("ADAPTIVE", "Adaptive", "Adaptive", "", 0),
        ("KELVIN", "Kelvin", "Kelvin", "", 1),
        ("CELSIUS", "Celsius", "Celsius", "", 2),
    ]
    temperature_unit_metric: EnumProperty(
        name='Temperature unit metric',
        items=temperature_units_metric,
        default='CELSIUS',
        get=get_temperature_unit_metric, set=set_temperature_unit_metric,
        #update=changeSubTypeFloat
    ) # type:  ignore

    def get_temperature_unit_imperial(self):
        temperature_unit = bpy.data.scenes['Scene'].unit_settings.temperature_unit
        res = [e for e in self.temperature_units_imperial if e[0]==temperature_unit][0][4]
        return res

    def set_temperature_unit_imperial(self, value):
        bpy.data.scenes['Scene'].unit_settings.temperature_unit = self.temperature_units_imperial[value][0]
        pass

    temperature_units_imperial = [
        ("ADAPTIVE", "Adaptive", "Adaptive", "", 0),
        ("KELVIN", "Kelvin", "Kelvin", "", 1),
        ("FAHRENHEIT", "Fahrenheit", "Fahrenheit", "", 2),
    ]
    temperature_unit_imperial: EnumProperty(
        name='Temperature unit imperial',
        items=temperature_units_imperial,
        default='FAHRENHEIT',
        get=get_temperature_unit_imperial, set=set_temperature_unit_imperial,
        #update=changeSubTypeFloat
    ) # type:  ignore

    def get_time_unit(self):
        time_unit = bpy.data.scenes['Scene'].unit_settings.time_unit
        res = [e for e in self.time_units if e[0]==time_unit][0][4]
        return res

    def set_time_unit(self, value):
        bpy.data.scenes['Scene'].unit_settings.time_unit = self.time_units[value][0]
        pass

    time_units = [
        ("ADAPTIVE", "Adaptive", "Adaptive", "", 0),
        ("DAYS", "Days", "Days", "", 1),
        ("HOURS", "Hours", "Hours", "", 2),
        ("MINUTES", "Minutes", "Minutes", "", 3),
        ("SECONDS", "Seconds", "Seconds", "", 4),
        ("MILLISECONDS", "Milliseconds", "Milliseconds", "", 5),
        ("MICROSECONDS", "Microseconds", "Microseconds", "", 6),
    ]
    time_unit: EnumProperty(
        name='Time unit',
        items=time_units,
        default='SECONDS',
        get=get_time_unit, set=set_time_unit,
        #update=changeSubTypeFloat
    ) # type:  ignore

    # color_pointer: bpy.props.PointerProperty(
    #     name="Object Reference",
    #     #poll=filter_kinds,  # seld.object_kinds can be "MESH" or "MESH,CURVE,.."
    #     #type=bpy.types.Object, # what kind of objects are we showing
    #     type=bpy.types.ColorSequence, # what kind of objects are we showing
    #     #update=process_from_socket,
    # ) # type: ignore

    # def quaternion_euler_order_update(self, context):
    #     #self.quaternion_euler_order = [e[0] for e in quaternion_euler_orders if e[3]==value][0]
    #     #self.quaternion_euler_order = value
    #     for I, elem in enumerate(self.quaternion_list_items):
    #         elem.quaternion_euler_order = self.quaternion_euler_order
    #         pass
    #     updateNode(self, context)
    #     pass

    quaternion_modes = [
        ("WXYZ", "WXYZ", "Convert components into quaternion", 0),
        ("SCALARVECTOR", "Scalar Vector", "Convert Scalar & Vector into quaternion", 1),
        ("EULER", "Euler Angles XYZ", "Convert Euler angles into quaternion", 2),
        ("AXISANGLE", "Angle Axis (Angle, XYZ)", "Convert Angle & Axis into quaternion", 3),
        #("MATRIX", "Matrix", "Convert Rotation Matrix into quaternion", 4),
    ]
    quaternion_mode : EnumProperty(
        name='Mode', description='The input component format of the quaternion',
        items=quaternion_modes,
        default="WXYZ",
        update=updateNode,
    ) # type: ignore

    quaternion_WXYZ_SCALARVECTOR_normalize : BoolProperty(
        name = "Normalize",
        description = "Normalize the output quaternion",
        default=False,
        update=updateNode,
    ) # type: ignore
    
    quaternion_euler_order: EnumProperty(
        name="Euler Order",
        description="Order of the Euler rotations",
        default="XYZ",
        items=quaternion_euler_orders,
        update=updateNode,
        #update=quaternion_euler_order_update
        #update=quaternion_euler_order_update
    ) # type: ignore

    matrix_modes1 = [
        ("NONE", "None", "None", 0),
        #("SCALARVECTOR", "Scalar Vector", "Convert Scalar & Vector into quaternion", 1),
        ("EULER", "Euler Angles XYZ", "Convert Euler angles into quaternion", 2),
        ("AXISANGLE", "Angle Axis (Angle, XYZ)", "Convert Angle & Axis into quaternion", 3),
    ]
    matrix_mode1 : EnumProperty(
        name='Matrix mode', description='The input component format of the matrix',
        items=matrix_modes1,
        default="NONE",
        update=updateNode,
    ) # type: ignore

    matrix_euler_order: EnumProperty(
        name="Euler Order",
        description="Order of the Euler rotations",
        default="XYZ",
        items=quaternion_euler_orders,
        update=updateNode,
    ) # type: ignore

    base_name = 'data '
    multi_socket_type = 'SvStringsSocket'

    use_alpha : BoolProperty(
        name = "Use Alpha",
        description = "Use Alpha channel in the Color",
        default=True,
        update=updateNode,
    ) # type: ignore


    def set_output_socketype(self, slot_bl_idnames):
        """
        1) if the input sockets are a mixed bag of bl_idnames we convert the output socket
        to a generic SvStringsSocket type
        2) if all input sockets where sv_get is successful are of identical bl_idname
        then set the output socket type to match that.
        3) no op if current output socket matches proposed new socket type.
        """

        if not slot_bl_idnames:
            return

        num_bl_idnames = len(set(slot_bl_idnames))
        new_socket_type = slot_bl_idnames[0] if num_bl_idnames == 1 else "SvStringsSocket"

        if self.outputs[0].bl_idname != new_socket_type:
            self.outputs[0].replace_socket(new_socket_type)

    def sv_init(self, context):
        self.width = 260

        self.inputs.new('SvStringsSocket', 'mask').label = "Mask"
        self.inputs['mask'].custom_draw = 'draw_mask_in_socket'

        if self.mode == 'BOOL_LIST_MODE':
            self.outputs.new('SvStringsSocket', 'data_output').label = 'Boolean'
        elif self.mode == 'INT_LIST_MODE':
            self.outputs.new('SvStringsSocket', 'data_output').label = 'Integers'
        elif self.mode == 'FLOAT_LIST_MODE':
            self.outputs.new('SvStringsSocket', 'data_output').label = 'Floats'
        elif self.mode == 'VECTOR_LIST_MODE':
            self.outputs.new('SvVerticesSocket', 'data_output').label = 'Vectors'
        elif self.mode == 'QUATERNION_LIST_MODE':
            self.outputs.new('SvQuaternionSocket', 'data_output').label = 'Quaternions'
        elif self.mode == 'MATRIX_LIST_MODE':
            self.outputs.new('SvMatrixSocket', 'data_output').label = 'Matrixes'
        elif self.mode == 'COLOR_LIST_MODE':
            self.outputs.new('SvColorSocket', 'data_output').label = 'Colors'
        elif self.mode == 'STRING_LIST_MODE':
            self.outputs.new('SvTextSocket', 'data_output').label = 'Strings'
        else:
            raise Exception(f"[func: sv_init] unknown mode {self.mode}.")
        Correct_ListInput_Length(self, context)
        
        pass

    def draw_buttons_ext(self, context, layout):
        r_unit_system = layout.row().split(factor=0.25)
        r_unit_system.column().label(text="Unit system:")
        r_unit_system.row().prop(self, "unit_system", expand=True)
        pass

    def draw_buttons(self, context, layout):
        cm_split = layout.row(align=True).split(factor=0.7, align=True)
        cm_split_c1 = cm_split.column(align=True)
        if self.mode == 'BOOL_LIST_MODE':
            cm_split_c1.prop(self, "bool_list_counter", text="List Length (bool)")
        elif self.mode == 'INT_LIST_MODE':
            cm_split_c1.prop(self, "int_list_counter", text="List Length (int)")
        elif self.mode == 'FLOAT_LIST_MODE':
            cm_split_c1.prop(self, "float_list_counter", text="List Length (float)")
        elif self.mode == 'VECTOR_LIST_MODE':
            cm_split_c1.prop(self, "vector_list_counter", text="List Length (vectors)")
        elif self.mode == 'QUATERNION_LIST_MODE':
            cm_split_c1.prop(self, "quaternion_list_counter", text="List Length (quaternions)")
        elif self.mode == 'MATRIX_LIST_MODE':
            cm_split_c1.prop(self, "matrix_list_counter", text="List Length (matrixes)")
        elif self.mode == 'COLOR_LIST_MODE':
            cm_split_c1.prop(self, "color_list_counter", text="List Length (colors)")
        elif self.mode == 'STRING_LIST_MODE':
            cm_split_c1.prop(self, "string_list_counter", text="List Length (strings)")
        else:
            raise Exception(f"[func: draw_buttons] unknown mode {self.mode}.")

        cm_split.column(align=True).prop(self, "mode", expand=False, text='')

        r_subtype_split1 = layout.row().split(factor=0.3)
        r_subtype_split1.column().label(text="Type:")
        r_subtype_split2 = r_subtype_split1.column().split(factor=0.6)
        if self.mode=='FLOAT_LIST_MODE':
            r_subtype_split2.column().prop(self, "subtype_float", expand=False, text='')
        elif self.mode=='VECTOR_LIST_MODE':
            r_subtype_split2.row().prop(self, "subtype_vector", expand=False, text='')
        elif self.mode=='COLOR_LIST_MODE':
            r_subtype_split2.row().prop(self, "subtype_color", expand=False, text='')
        elif self.mode=='STRING_LIST_MODE':
            r_subtype_split2.row().prop(self, "subtype_string", expand=False, text='')
        elif self.mode=='QUATERNION_LIST_MODE':
            label_row = r_subtype_split2.row()
            label_row.label(text='Quaternion')
        elif self.mode=='MATRIX_LIST_MODE':
            if self.matrix_mode1=='NONE':
                pass
            else:
                label_row = r_subtype_split2.row()
                #label_row.label(text='Matrix')
                if self.unit_system=='METRIC':
                    label_row.prop(self, 'length_unit_metric', text='')
                elif self.unit_system=='IMPERIAL':
                    label_row.prop(self, 'length_unit_imperial', text='')
        else:
            r_subtype_split1.enabled = False
            pass

        if self.unit_system=='METRIC':
            if  self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='DISTANCE' or \
                self.mode=='VECTOR_LIST_MODE' and (self.subtype_vector=='TRANSLATION' or self.subtype_vector=='XYZ_LENGTH'):
                    r_subtype_split2.column().prop(self, 'length_unit_metric', expand=False, text='')
            elif self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='ANGLE' or \
                self.mode=='VECTOR_LIST_MODE' and (self.subtype_vector=='EULER' or self.subtype_vector=='AXISANGLE' or self.subtype_vector=='DIRECTION') or \
                self.mode=='QUATERNION_LIST_MODE' and (self.quaternion_mode=='EULER' or self.quaternion_mode=='AXISANGLE') or \
                self.mode=='MATRIX_LIST_MODE' and (self.matrix_mode1=='EULER' or self.matrix_mode1=='AXISANGLE'):
                    r_subtype_split2.column().prop(self, 'system_rotation', expand=False, text='')
            elif self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='TEMPERATURE':
                    r_subtype_split2.column().prop(self, 'temperature_unit_metric', expand=False, text='')
            elif self.mode=='COLOR_LIST_MODE' and (self.subtype_color=='COLOR' or self.subtype_color=='COLOR_GAMMA'):
                    r_subtype_split2.column().prop(self, 'use_alpha')
            else:
                    r_subtype_split2.column().label(text='')
        elif self.unit_system=='IMPERIAL':
            if  self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='DISTANCE' or \
                self.mode=='VECTOR_LIST_MODE' and (self.subtype_vector=='TRANSLATION' or self.subtype_vector=='XYZ_LENGTH'):
                    r_subtype_split2.column().prop(self, 'length_unit_imperial', expand=False, text='')
            elif self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='ANGLE' or \
                 self.mode=='VECTOR_LIST_MODE' and (self.subtype_vector=='EULER' or self.subtype_vector=='AXISANGLE' or self.subtype_vector=='DIRECTION') or \
                 self.mode=='QUATERNION_LIST_MODE' and (self.quaternion_mode=='EULER' or self.quaternion_mode=='AXISANGLE') or \
                 self.mode=='MATRIX_LIST_MODE' and (self.matrix_mode1=='EULER' or self.matrix_mode1=='AXISANGLE'):
                    r_subtype_split2.column().prop(self, 'system_rotation', expand=False, text='')
            elif self.mode=='FLOAT_LIST_MODE' and self.subtype_float=='TEMPERATURE':
                    r_subtype_split2.column().prop(self, 'temperature_unit_imperial', expand=False, text='')
            elif self.mode=='COLOR_LIST_MODE' and (self.subtype_color=='COLOR' or self.subtype_color=='COLOR_GAMMA'):
                    r_subtype_split2.column().prop(self, 'use_alpha')
            else:
                    r_subtype_split2.column().label(text='')
        else:
            if   self.mode=='FLOAT_LIST_MODE' and (self.subtype_float=='ANGLE') or \
                 self.mode=='VECTOR_LIST_MODE' and (self.subtype_vector=='EULER' or self.subtype_vector=='AXISANGLE' or self.subtype_vector=='DIRECTION') or \
                 self.mode=='QUATERNION_LIST_MODE' and (self.quaternion_mode=='EULER' or self.quaternion_mode=='AXISANGLE') or \
                 self.mode=='MATRIX_LIST_MODE' and (self.matrix_mode1=='EULER' or self.matrix_mode1=='AXISANGLE'):
                    r_subtype_split2.column().prop(self, 'system_rotation', expand=False, text='')
            elif self.mode=='COLOR_LIST_MODE' and (self.subtype_color=='COLOR' or self.subtype_color=='COLOR_GAMMA'):
                    r_subtype_split2.column().prop(self, 'use_alpha')
            else:
                    r_subtype_split2.column().label(text='')

        if self.mode=='QUATERNION_LIST_MODE':
            qp_row = layout.row() # quaternion_params_row
            qp_row_s = qp_row.split(factor=0.3)
            qp_row_s_c1 = qp_row_s.column()
            qp_row_s_c1.label(text='Quaternion mode:')
            qp_row_s_s = qp_row_s.column().split(factor=0.6)
            qp_row_s_s_c1 = qp_row_s_s.column()
            qp_row_s_s_c1.alignment="LEFT"
            qp_row_s_s_c1.prop(self, 'quaternion_mode', text='')
            qp_row_s_s_c2 = qp_row_s_s.row()
            if self.quaternion_mode=='WXYZ' or self.quaternion_mode=='SCALARVECTOR':
                qp_row_s_s_c2.prop(self, 'quaternion_WXYZ_SCALARVECTOR_normalize', toggle=True)
            elif self.quaternion_mode=='EULER':
                qp_row_s_s_c2.prop(self, 'quaternion_euler_order', text='')
            else:
                qp_row_s_s_c2.label(text='')
            pass
        elif self.mode=='MATRIX_LIST_MODE':
            # Angles
            qp_row = layout.row() # matrix_params_row
            qp_row_s = qp_row.split(factor=0.3)
            qp_row_s_c1 = qp_row_s.column()
            qp_row_s_c1.label(text='Matrix mode:')
            qp_row_s_s = qp_row_s.column().split(factor=0.6)
            qp_row_s_s_c1 = qp_row_s_s.column()
            qp_row_s_s_c1.alignment="LEFT"
            qp_row_s_s_c1.prop(self, 'matrix_mode1', text='')
            qp_row_s_s_c2 = qp_row_s_s.row()
            if self.matrix_mode1=='NONE':
                pass
            elif self.matrix_mode1=='EULER':
                qp_row_s_s_c2.prop(self, 'matrix_euler_order', text='')
            else:
                qp_row_s_s_c2.label(text='')
            pass

        # invert_mask_row = layout.row()
        # invert_mask_row.prop(self, "invert_mask", text="Invert mask", toggle=True)
        grid_service_operators = layout.grid_flow(row_major=True, columns=10, align=True)
        
        if self.mode=='STRING_LIST_MODE':
            self.wrapper_tracked_ui_draw_op(grid_service_operators, SvUpdateTextInListInputNode.bl_idname, text='update text', icon='FILE_REFRESH')
        else:
            grid_service_operators.label(text='')
        invert_mask_prop = grid_service_operators.column()
        # if mask socked is connected then do not show source list controls
        if self.inputs["mask"].is_linked==True:
            invert_mask_prop.enabled = False
            invert_mask_prop.prop(self, "invert_mask", text='', icon='UV_SYNC_SELECT', emboss = False)
        else:
            invert_mask_prop.prop(self, "invert_mask", text='', icon='UV_SYNC_SELECT')

        #layout.row().prop(self, "copy_clipboard", text="Copy data to clipboard", toggle=True)
        self.wrapper_tracked_ui_draw_op(grid_service_operators, SvCopyTextToClipboard.bl_idname, text='', icon='COPYDOWN')
        self.wrapper_tracked_ui_draw_op(grid_service_operators, SvPasteTextFromClipboard.bl_idname, text='', icon='PASTEDOWN')
        pass

       
        align = True
        if self.mode=='VECTOR_LIST_MODE' or self.mode=='COLOR_LIST_MODE':
            align=False

        col = layout.column(align=align)
        J=0
        if self.mode == 'BOOL_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=4, align=True, even_rows=False)
            for I, elem in enumerate(self.bool_list_items):
                c1 = grid_row.column()
                c2 = grid_row.column()
                c3 = grid_row.column()
                c4 = grid_row.column()

                c1.ui_units_x = 50
                c2.ui_units_x = 5
                c3.ui_units_x = 15
                c4.ui_units_x = 30
                
                c4.alignment = 'RIGHT'

                c1.prop(elem, f'elem', text=str(I))
                # if mask socked is connected then do not show source list controls
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                        value_label = ''
                    else:
                        index_label = f'{J}'
                        value_label = f'{elem.elem}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                else:
                    index_label = f'{J}'
                    value_label = ''
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                pass

        elif self.mode == 'INT_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=4, align=True, even_rows=False)
            for I, elem in enumerate(self.int_list_items):
                c1 = grid_row.column()
                c2 = grid_row.column()
                c3 = grid_row.column()
                c4 = grid_row.column()

                c1.ui_units_x = 50
                c2.ui_units_x = 5
                c3.ui_units_x = 15
                c4.ui_units_x = 30
                
                c4.alignment = 'RIGHT'

                c1.prop(elem, f'elem', text=str(I))
                # if mask socked is connected then do not show source list controls
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                        value_label = ''
                    else:
                        index_label = f'{J}'
                        value_label = f'{elem.elem}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                else:
                    index_label = f'{J}'
                    value_label=''
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                pass

        elif self.mode == 'FLOAT_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=4, align=True, even_rows=False)
            for I, elem in enumerate(self.float_list_items):
                c1 = grid_row.column()
                c2 = grid_row.column()
                c3 = grid_row.column()
                c4 = grid_row.column()

                c1.ui_units_x = 50
                c2.ui_units_x = 5
                c3.ui_units_x = 15
                c4.ui_units_x = 30

                c4.alignment = 'RIGHT'

                c1.prop(elem, self.subtype_float, text=str(I))

                # if mask socked is connected then do not show source list controls
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                        value_label = ''
                    else:
                        index_label = f'{J}'
                        value_label = f'{elem.elem}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                else:
                    index_label = f'{J}'
                    value_label = ''
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                    c4.label(text=value_label)
                pass

        elif self.mode == 'VECTOR_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=3, align=True, even_rows=False)
            for I, elem in enumerate(self.vector_list_items):
                c1 = grid_row.column()
                c2 = grid_row.column()
                c3 = grid_row.column()

                c1.ui_units_x = 85
                c2.ui_units_x = 5
                c3.ui_units_x = 10

                if self.subtype_vector=='DIRECTION':
                    c1_r = c1.row(align=True)
                    c1_r.column(align=True).prop(elem, self.subtype_vector, icon_only=True)
                    c1_r.column(align=True).prop(elem, 'EULER', icon_only=True)
                    pass
                else:
                    c1.row().prop(elem, self.subtype_vector, icon_only=True)

                # if mask socked is connected then do not show source list controls
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                    else:
                        index_label = f'{J}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                else:
                    index_label = f'{J}'
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                pass
            
        elif self.mode == 'QUATERNION_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=3, align=True, even_rows=False)
            for I, elem in enumerate(self.quaternion_list_items):
                c1 = grid_row.row()
                c2 = grid_row.column()
                c3 = grid_row.column()

                c1.ui_units_x = 85
                c2.ui_units_x = 5
                c3.ui_units_x = 10

                c1.prop(elem, self.quaternion_mode+'_UI', icon_only=True)
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                    else:
                        index_label = f'{J}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                else:
                    index_label = f'{J}'
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                pass
            
        elif self.mode == 'MATRIX_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=3, align=True, even_rows=False)
            for I, elem in enumerate(self.matrix_list_items):
                c1 = grid_row.row()
                c2 = grid_row.column()
                c3 = grid_row.column()

                c1.ui_units_x = 85
                c2.ui_units_x = 5
                c3.ui_units_x = 10

                c1_r = c1.row()
                if self.matrix_mode1=='NONE':
                    matrix_grid = c1_r.grid_flow(row_major=True, columns=4, align=True)
                    # for i in range(4):
                    #     for j in range(i, 16, 4):
                    #         matrix_grid.prop(elem, 'MATRIX_UI', text='', index=j)
                    #     pass
                    for I in range(16):
                        matrix_grid.prop(elem, 'MATRIX_UI', text='', index=I)
                    pass
                    # for i in range(16):
                    #     matrix_grid.prop(elem, 'elem', text='', index=i)
                    # pass
                else:
                    if self.matrix_mode1=='EULER':
                        matrix_grid = c1_r.grid_flow(row_major=False, columns=3, align=True)
                        matrix_grid.label(text='Location:')
                        if self.unit_system=='NONE':
                            matrix_grid.prop(elem, 'EULER_LOCATION_UI_NONE', text='')
                        else:
                            matrix_grid.prop(elem, 'EULER_LOCATION_UI_TRANSLATION', text='')

                        matrix_grid.label(text='Scale:')
                        matrix_grid.prop(elem, 'EULER_SCALE_UI', text='')
                        
                        matrix_grid.label(text='Angle:')
                        matrix_grid.prop(elem, 'EULER_ANGLE_UI', text='')
                        pass
                    elif self.matrix_mode1=='AXISANGLE':
                        matrix_grid = c1_r.grid_flow(row_major=False, columns=3, align=True)
                        matrix_grid.label(text='Location:')
                        if self.unit_system=='NONE':
                            matrix_grid.prop(elem, 'AXISANGLE_LOCATION_UI_NONE', text='')
                        else:
                            matrix_grid.prop(elem, 'AXISANGLE_LOCATION_UI_TRANSLATION', text='')

                        matrix_grid.label(text='Scale:')
                        matrix_grid.prop(elem, 'AXISANGLE_SCALE_UI', text='')
                        
                        matrix_grid.label(text='Axis:')
                        matrix_grid.prop(elem, 'AXISANGLE_VECTOR_UI', text='')

                        c1_r_r = grid_row.row(align=True)
                        c1_r_r.column().label(text='Angle:')
                        c1_r_r.column().prop(elem, 'AXISANGLE_ANGLE_UI', text='')
                        grid_row.column()
                        grid_row.column()
                        pass
                    else:
                        raise Exception(f'Unknown Matrix mode: {self.matrix_mode1}')
                
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                    else:
                        index_label = f'{J}'
                        J+=1
                    if self.matrix_mode1=='NONE':
                        pass
                    else:
                        c2.row().label(text='')
                        c3.row().label(text='')
                    c2.row().prop(elem, f'item_enable', icon_only=True)
                    c3.row().label(text=index_label)
                else:
                    index_label = f'{J}'
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)

                col.row()
                pass

        elif self.mode == 'COLOR_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=3, align=True, even_rows=False)
            for I, elem in enumerate(self.color_list_items):
                c1 = grid_row.column().row(align=True) # row.row - to align=True only one row. If only row() then all elems in grid will connect each other.
                c2 = grid_row.column()
                c3 = grid_row.column()

                c1.ui_units_x = 85
                c2.ui_units_x = 5
                c3.ui_units_x = 10

                if self.subtype_color=='COLOR' or self.subtype_color=='COLOR_GAMMA':
                    if self.use_alpha==True:
                        c1.row(align=True).prop(elem, self.subtype_color, icon_only=True)
                        c1.row(align=True).prop(elem, 'elem', icon_only=True)
                    else:
                        c1.row(align=True).prop(elem, self.subtype_color+'_WITHOUT_ALPHA', icon_only=True)
                        c1.row(align=True).prop(elem, 'elem_WITHOUT_ALPHA', icon_only=True)
                    
                    pass
                else:
                    c1.prop(elem, self.subtype_color, icon_only=True)

                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                    else:
                        index_label = f'{J}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                else:
                    index_label = f'{J}'
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                pass
        elif self.mode == 'STRING_LIST_MODE':
            grid_row = col.row().grid_flow(row_major=True, columns=3, align=True, even_rows=False)
            for I, elem in enumerate(self.string_list_items):

                c1 = grid_row.column()
                c2 = grid_row.column()
                c3 = grid_row.column()

                c1.ui_units_x = 85
                c2.ui_units_x = 5
                c3.ui_units_x = 10

                c1.prop(elem, self.subtype_string, icon_only=False, text='')
                if self.inputs["mask"].is_linked==False:
                    if elem.item_enable==False:
                        index_label = '-'
                    else:
                        index_label = f'{J}'
                        J+=1
                    c2.prop(elem, f'item_enable', icon_only=True)
                    c3.label(text=index_label)
                else:
                    index_label = f'{J}'
                    J+=1
                    c2.label(text='')
                    c3.label(text=index_label)
                pass
        else:
            raise Exception(f"[func: draw_buttons] unknown mode {self.mode}.")

    def draw_buttons_3dpanel(self, layout, in_menu=None):
        if not in_menu:
            menu = layout.row(align=True).operator('node.popup_3d_menu', text=f'Show: "{self.label or self.name}"')
            menu.tree_name = self.id_data.name
            menu.node_name = self.name
        else:
            layout.label(text=self.label or self.name)
            if self.mode == 'vector':
                colum_list = layout.column(align=True)
                for i in range(self.v_int):
                    row = colum_list.row(align=True)
                    for j in range(3):
                        row.prop(self, 'vector_list', index=i*3+j, text='XYZ'[j]+(self.label if self.label else self.name))
            else:
                colum_list = layout.column(align=True)
                for i in range(self.int_list_counter):
                    row = colum_list.row(align=True)
                    row.prop(self, self.mode, index=i, text=str(i)+(self.label if self.label else self.name))
                    row.scale_x = 0.8

    def updateTextData(self, context):
        self.process_node(context)
        pass

    def process(self):
        mask_in_socket = self.inputs['mask'] #.sv_get(deepcopy=False)
        if mask_in_socket.is_linked==False:
            mask_in_data = [[]]
        else:
            mask_in_data = mask_in_socket.sv_get(deepcopy=False)
        mask_in_data = ensure_nesting_level(mask_in_data, 2)

        data = []
        append = data.append
        if self.outputs[0].is_linked:
            if mask_in_socket.is_linked==False:

                if self.mode == 'BOOL_LIST_MODE':
                    lst = [[elem.elem for elem in self.bool_list_items if elem.item_enable==True]]

                elif self.mode == 'INT_LIST_MODE':
                    lst = [[elem.elem for elem in self.int_list_items if elem.item_enable==True]]

                elif self.mode == 'FLOAT_LIST_MODE':
                    lst = [[elem.elem for elem in self.float_list_items if elem.item_enable==True]]

                elif self.mode == 'VECTOR_LIST_MODE':
                    lst = [[tuple(Vector(elem.elem)) for elem in self.vector_list_items if elem.item_enable==True]]

                elif self.mode == 'QUATERNION_LIST_MODE':
                    if (self.quaternion_mode=='WXYZ' or self.quaternion_mode=='SCALARVECTOR') and self.quaternion_WXYZ_SCALARVECTOR_normalize==True:
                        lst = [[Quaternion(elem.elem).normalized() for elem in self.quaternion_list_items if elem.item_enable==True]]
                    elif self.quaternion_mode=='EULER':
                        lst = [[Quaternion(elem.elem).to_euler(self.quaternion_euler_order).to_quaternion() for elem in self.quaternion_list_items if elem.item_enable==True]]
                    else:
                        lst = [[Quaternion(elem.elem) for elem in self.quaternion_list_items if elem.item_enable==True]]

                elif self.mode == 'MATRIX_LIST_MODE':
                    #lst = [[Matrix(elem.MATRIX_UI) for elem in self.matrix_list_items if elem.item_enable==True]]
                    lst = [elem.MATRIX[:] for elem in self.matrix_list_items if elem.item_enable==True]
                    lst = [[Matrix(m) for m in np.transpose( np.array(list(chain(*lst),)).reshape(-1,4,4), (0,1,2) ).tolist()]]

                elif self.mode == 'COLOR_LIST_MODE':
                    if self.use_alpha==True:
                        lst = [[tuple(elem.elem) for elem in self.color_list_items if elem.item_enable==True]]
                    else:
                        lst = [[tuple(elem.elem[:3]) for elem in self.color_list_items if elem.item_enable==True]]

                elif self.mode == 'STRING_LIST_MODE':
                    lst = [[elem.elem for elem in self.string_list_items if elem.item_enable==True]]
                
                else:
                    raise Exception(f"[func: process] unknown mode {self.mode}.")
                
                data = lst

            else:

                if self.mode == 'BOOL_LIST_MODE':
                    lst = [[elem.elem for elem in self.bool_list_items]]
                elif self.mode == 'INT_LIST_MODE':
                    lst = [[elem.elem for elem in self.int_list_items]]
                elif self.mode == 'FLOAT_LIST_MODE':
                    lst = [[elem.elem for elem in self.float_list_items]]
                elif self.mode == 'VECTOR_LIST_MODE':
                    lst = [[tuple(Vector(elem.elem)) for elem in self.vector_list_items]]
                elif self.mode == 'QUATERNION_LIST_MODE':
                    if (self.quaternion_mode=='WXYZ' or self.quaternion_mode=='SCALARVECTOR') and self.quaternion_WXYZ_SCALARVECTOR_normalize==True:
                        lst = [[Quaternion(elem.elem).normalized() for elem in self.quaternion_list_items]]
                    elif self.quaternion_mode=='EULER':
                        lst = [[Quaternion(elem.elem).to_euler(self.quaternion_euler_order).to_quaternion() for elem in self.quaternion_list_items]]
                    else:
                        lst = [[Quaternion(elem.elem) for elem in self.quaternion_list_items]]
                elif self.mode == 'MATRIX_LIST_MODE':
                    lst = [[Matrix(elem.MATRIX) for elem in self.matrix_list_items]]
                elif self.mode == 'COLOR_LIST_MODE':
                    if self.use_alpha==True:
                        lst = [[tuple(elem.elem) for elem in self.color_list_items]]
                    else:
                        lst = [[tuple(elem.elem[:3]) for elem in self.color_list_items]]
                elif self.mode == 'STRING_LIST_MODE':
                    lst = [[elem.elem for elem in self.string_list_items]]
                else:
                    raise Exception(f"[func: process] unknown mode {self.mode}.")
                
                lst_0 = lst[0]
                for mask in mask_in_data:
                    # if mask is zero or not connected then do not mask any. Except of inversion,
                    if not mask:
                        np_mask = np.zeros(len(lst_0), dtype=bool)
                        if self.mask_inversion==True:
                            np_mask = np.invert(np_mask)
                        mask = np_mask.tolist()
                    else:
                        if self.mask_mode=='BOOLEAN':
                            if self.mask_inversion==True:
                                mask = list( map( lambda v: False if v==0 else True, mask) )
                                mask = mask[:len(lst_0)]
                                np_mask = np.zeros(len(lst_0), dtype=bool)
                                np_mask[0:len(mask)]=mask
                                np_mask = np.invert(np_mask)
                                mask = np_mask.tolist()
                            pass
                        elif self.mask_mode=='INDEXES':
                            mask_len = len(lst_0)
                            mask_range = []
                            for x in mask:
                                if type(x)==bool:
                                    if x==False:
                                        x=0
                                    else:
                                        x=1
                                if -mask_len<x<mask_len:
                                    mask_range.append(x)
                            np_mask = np.zeros(len(lst_0), dtype=bool)
                            np_mask[mask_range] = True
                            if self.mask_inversion==True:
                                np_mask = np.invert(np_mask)
                            mask = np_mask.tolist()
                    mask = [not(e) for e in mask]
                    lst_masked = list(compress(lst_0, mask))
                    data.append(lst_masked)
                    pass
        self.outputs[0].sv_set(data)

    def dataAsString(self):

        if self.mode == 'BOOL_LIST_MODE':
            lst = [[str(elem.elem) for elem in self.bool_list_items]]
        elif self.mode == 'INT_LIST_MODE':
            lst = [[str(elem.elem) for elem in self.int_list_items]]
        elif self.mode == 'FLOAT_LIST_MODE':
            lst = [[str(elem.elem) for elem in self.float_list_items]]
        elif self.mode == 'VECTOR_LIST_MODE':
            lst = [[str(tuple(Vector(elem.elem)))[1:-1] for elem in self.vector_list_items]]
        elif self.mode == 'MATRIX_LIST_MODE':
            lst = [[str(tuple(elem.elem))[1:-1] for elem in self.matrix_list_items]]
        elif self.mode == 'QUATERNION_LIST_MODE':
            lst = [[str(getattr(elem, f'{self.quaternion_mode}_UI')[:])[1:-1] for elem in self.quaternion_list_items]]
        elif self.mode == 'COLOR_LIST_MODE':
            if self.use_alpha==True:
                lst = [[str(tuple(elem.elem))[1:-1] for elem in self.color_list_items]]
            else:
                lst = [[str(tuple(elem.elem[:3]))[1:-1] for elem in self.color_list_items]]
        elif self.mode == 'STRING_LIST_MODE':
            lst = [[elem.elem for elem in self.string_list_items]]
        else:
            ShowMessageBox(f"In copy to the clipboard unknown mode: '{self.mode}'. Data not copied.", icon='ERROR')
            #raise Exception(f"[func: process] unknown mode {self.mode}.")
        
        lst_0 = lst[0]
        res = '\n'.join( [str(s) for s in lst_0] )

        return res
    
    def convertTextToData(self, text):
        text = text.strip()
        if self.mode == 'BOOL_LIST_MODE':
            text = text.replace("  ", " ")
            txt1 = re.sub(r'([^\,]) ([^$])', r'\1, \2', ' '.join(text.splitlines()))
            txt2 = re.sub(r' ', r'', txt1)
            txt2_arr = txt2.split(',')
            arr = []
            for elem in txt2_arr:
                if   elem in [1, '1', 'ON', 'on', 'On', 'TRUE', 'True', 'true', 'YES','Yes', 'yes',]:
                    c_elem = True
                elif elem in [0, '0', 'OFF', 'Off', 'off', 'FALSE', 'False', 'false', 'NO','No', 'no',]:
                    c_elem = False
                else:
                    c_elem = bool(elem)
                arr.append(c_elem)
                pass
            self.bool_list_counter = len(arr)
            for I, elem in enumerate(arr):
                self.bool_list_items[I].elem = elem
            pass

        elif self.mode == 'INT_LIST_MODE':
            # Convert multiline formats of vectors of floats as single line of floats
            text = text.replace("  ", " ")
            txt1 = re.sub(r'([^\,]) ', r'\1, ', ' '.join(text.splitlines()))
            txt2 = re.sub(r' ', r'', txt1)
            txt2_arr = txt2.split(',')
            arr = []
            for I, elem in enumerate(txt2_arr):
                try:
                    c_elem = int(float(elem))
                except ValueError:
                    ShowMessageBox(f"In paste the int numbers raise error: value '{str(elem)}' is not an integer number. Data not pasted.", icon='ERROR')
                    return
                arr.append(c_elem)
                pass
            self.int_list_counter = len(arr)
            for I, elem in enumerate(arr):
                self.int_list_items[I].elem = elem
            pass

        elif self.mode == 'FLOAT_LIST_MODE':
            # Convert multiline formats of vectors of floats as single line of floats
            text = text.replace("  ", " ")
            txt1 = re.sub(r'([^\,]) ', r'\1, ', ' '.join(text.splitlines()))
            txt2 = re.sub(r' ', r'', txt1)
            txt2_arr = txt2.split(',')
            arr = []
            for I, elem in enumerate(txt2_arr):
                try:
                    c_elem = float(elem)
                except ValueError:
                    ShowMessageBox(f"In paste the float numbers raise error: value '{str(elem)}' is not an integer number. Data not pasted.", icon='ERROR')
                    return
                arr.append(c_elem)
                pass
            self.float_list_counter = len(arr)
            for I, elem in enumerate(arr):
                self.float_list_items[I].elem = elem
            pass

        elif self.mode == 'VECTOR_LIST_MODE':
            text = text.replace("  ", " ")
            _arr_lines = text.splitlines()
            arr_lines = []
            for line in _arr_lines:
                line = line.strip()
                if line:
                    arr_lines.append(line)
                pass
            pass

            
            arr = []
            for line in arr_lines:
                txt1 = re.sub(r'([^\,]) ', r'\1, ', line)
                txt2 = re.sub(r' ', r'', txt1)
                txt2_arr = txt2.split(',')
                arr1 = []
                for I, elem in enumerate(txt2_arr):
                    try:
                        c_elem = float(elem)
                    except ValueError:
                        ShowMessageBox(f"In paste the float numbers raise error: value '{str(elem)}' is not an float number. Data not pasted.", icon='ERROR')
                        return
                    arr1.append(c_elem)
                    pass
                arr.append(arr1)
                pass
            
            self.vector_list_counter = len(arr)
            for I, line in enumerate(arr):
                self.vector_list_items[I].elem = [0]*len(self.vector_list_items[I].elem)
                for J, elem in enumerate(line):
                    if len(self.vector_list_items[I].elem)<=J:
                        break
                    self.vector_list_items[I].elem[J] = elem
            pass

        elif self.mode == 'COLOR_LIST_MODE':
            text = text.replace("  ", " ")
            _arr_lines = text.splitlines()
            arr_lines = []
            for line in _arr_lines:
                line = line.strip()
                if line:
                    arr_lines.append(line)
                pass
            pass

            
            arr = []
            for line in arr_lines:
                txt1 = re.sub(r'([^\,]) ', r'\1, ', line)
                txt2 = re.sub(r' ', r'', txt1)
                txt2_arr = txt2.split(',')
                arr1 = []
                for I, elem in enumerate(txt2_arr):
                    try:
                        c_elem = float(elem)
                    except ValueError:
                        ShowMessageBox(f"In paste the float numbers raise error: value '{str(elem)}' is not an float number. Data not pasted.", icon='ERROR')
                        return
                    arr1.append(c_elem)
                    pass
                arr.append(arr1)
                pass
            
            self.color_list_counter = len(arr)
            for I, line in enumerate(arr):
                self.color_list_items[I].elem = [0]*len(self.color_list_items[I].elem)
                self.color_list_items[I].elem[-1]=1.0
                for J, elem in enumerate(line):
                    if len(self.color_list_items[I].elem)<=J:
                        break
                    self.color_list_items[I].elem[J] = elem
            pass

        elif self.mode == 'QUATERNION_LIST_MODE':
            text = text.replace("  ", " ")
            _arr_lines = text.splitlines()
            arr_lines = []
            for line in _arr_lines:
                line = line.strip()
                if line:
                    arr_lines.append(line)
                pass
            pass

            
            arr = []
            for line in arr_lines:
                txt1 = re.sub(r'([^\,]) ', r'\1, ', line)
                txt2 = re.sub(r' ', r'', txt1)
                txt2_arr = txt2.split(',')
                arr1 = []
                for I, elem in enumerate(txt2_arr):
                    try:
                        c_elem = float(elem)
                    except ValueError:
                        ShowMessageBox(f"In paste the float numbers raise error: value '{str(elem)}' is not an float number. Data not pasted.", icon='ERROR')
                        return
                    arr1.append(c_elem)
                    pass
                arr.append(arr1)
                pass
            
            self.quaternion_list_counter = len(arr)
            for I, line in enumerate(arr):
                attr = getattr(self.quaternion_list_items[I], f'{self.quaternion_mode}_UI')
                attr[:] = [0]*len(attr)
                for J, elem in enumerate(line):
                    if len(attr)<=J:
                        break
                    attr[J] = elem
            pass

        elif self.mode == 'STRING_LIST_MODE':
            text = text.replace("  ", " ")
            _arr_lines = text.splitlines()
            arr_lines = []
            for line in _arr_lines:
                line = line.strip()
                if line:
                    arr_lines.append(line)
                pass
            pass

            self.string_list_counter = len(arr_lines)
            for I, line in enumerate(arr_lines):
                self.string_list_items[I].elem = line
            pass

        else:
            mode_name = [m for m in self.modes if m[0]==self.mode][0][1]
            ShowMessageBox(f"Data NOT pasted from the Clipboard, for: '{mode_name}'", icon='ERROR')
            return

        ShowMessageBox(f"Data pasted from the Clipboard")
        pass


classes = [SvCopyTextToClipboard, SvPasteTextFromClipboard, SvUpdateTextInListInputNode, SvListInputBoolEntry, SvListInputIntEntry, SvListInputFloatEntry, SvListInputVectorEntry, SvListInputMatrixEntry, SvListInputQuaternionEntry, SvListInputColorEntry, SvListInputStringEntry, SvListInputNodeMK2]
register, unregister = bpy.utils.register_classes_factory(classes)
