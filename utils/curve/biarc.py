# This file is part of project Sverchok. It's copyrighted by the contributors
# recorded in the version control history of the file, available from
# its original location https://github.com/nortikin/sverchok/commit/master
#
# SPDX-License-Identifier: GPL3
# License-Filename: LICENSE

import numpy as np
from math import sin, cos, acos
from cmath import exp as cexp, phase as arg

from mathutils import Vector, Matrix

from sverchok.utils.curve.core import SvCurve
from sverchok.utils.curve.primitives import SvCircle
from sverchok.utils.curve.algorithms import SvConcatCurve

class SvBiArc(SvConcatCurve):
    def __init__(self, arc1, arc2):
        super().__init__([arc1, arc2])
        self.arc1 = arc1
        self.arc2 = arc2

    @staticmethod
    def calc(point_a, point_b, tangent_a, tangent_b, p):
        xx = point_b - point_a
        xx /= np.linalg.norm(xx)
        zz = np.cross(tangent_a, tangent_b)
        volume = np.dot(xx, zz)
        if abs(volume) > 1e-6:
            raise Exception(f"Provided tangents are not coplanar, volume={volume}, zz={zz}")
        c = np.linalg.norm(point_b - point_a) * 0.5
        origin = (point_a + point_b) * 0.5

        zz = np.cross(tangent_a, xx)
        
        tangent_a /= np.linalg.norm(tangent_a)
        tangent_b /= np.linalg.norm(tangent_b)
        zz /= np.linalg.norm(zz)
        yy = np.cross(zz, xx)

        norm_a = np.cross(zz, tangent_a)
        norm_b = - np.cross(zz, tangent_b)

        matrix_inv = np.stack((xx, yy, zz))
        matrix = np.linalg.inv(matrix_inv)

        # TODO: sin(arccos(...)) = ...
        alpha = acos(np.dot(tangent_a, xx))
        beta = acos(np.dot(tangent_b, xx))
        if np.dot(tangent_a, yy) > 0:
            alpha = - alpha
        if np.dot(tangent_b, yy) > 0:
            beta = - beta
        #print("A", alpha, "B", beta)

        omega = (alpha + beta) * 0.5
        r1 = c/(sin(alpha) + sin(omega) / p)
        r2 = c/(sin(beta) + p * sin(omega))
        #print("R1", r1, "R2", r2)
        theta1 = 2 * arg(cexp(-1j * alpha) + cexp(-1j*omega)/p)
        theta2 = 2 * arg(cexp(1j*beta) + p*cexp(1j*omega))
        #print("T1", theta1, "T2", theta2)

        vectorx_a = r1 * norm_a
        center1 = point_a + vectorx_a
        center2 = point_b + r2 * norm_b

        arc1 = SvCircle(#radius = r1,
                    center = center1,
                    normal = -zz, vectorx = point_a - center1)
        if theta1 > 0:
            arc1.u_bounds = (0.0, theta1)
        else:
            theta1 = -theta1
            arc1.u_bounds = (0.0, theta1)
            arc1.normal = -arc1.normal

        junction = arc1.evaluate(theta1)
        #print("J", junction)

        arc2 = SvCircle(#radius = r2,
                    center = center2,
                    normal = -zz, vectorx = junction - center2)
        if theta2 > 0:
            arc2.u_bounds = (0.0, theta2)
        else:
            theta2 = -theta2
            arc2.u_bounds = (0.0, theta2)
            arc2.normal = -arc2.normal

        return SvBiArc(arc1, arc2)

