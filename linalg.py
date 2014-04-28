import numpy as np
import operator
import numbers
import math

def rotation_vector(v1, v2):
    cross = v1 ** v2

    # If the vectors are pointing in opposite
    # directions, we will get the zero vector,
    # which is not useful. We just need to find
    # some vector which is perpendicular to either
    # vector. There are an infinite number of them,
    # but it doesn't matter which one we pick.
    if all(a == 0 for a in cross):
        basis = np.array([0, 0, 0])
        min_index = min(enumerate(v1), key=lambda x: x[1])[0]
        basis[min_index] = 1
        cross = v1 ** basis

    return cross.normalize()

def _make_op(op):
    def func(self, other):
        if isinstance(other, Vector):
            other = other.v

        result = op(self.v, other)
        try:
            return Vector(*result)
        except TypeError:
            return np.asscalar(result)

    return func

def _make_rop(op):
    return _make_op(lambda a, b: op(b, a))

class Vector:
    symbols = ['x', 'y', 'z']
    sym_to_index = {sym: i for i, sym in enumerate(symbols)}

    def __init__(self, *vals):
        self.v = np.array(vals)

    __mul__ = _make_op(np.dot)
    __rmul__ = _make_rop(np.dot)

    __pow__ = _make_op(np.cross)
    __rpow__ = _make_rop(np.cross)

    __add__ = _make_op(operator.add)
    __radd__ = __add__

    __sub__ = _make_op(operator.sub)
    __rsub__ = _make_rop(operator.sub)

    __truediv__ = _make_op(operator.truediv)
    __rtruediv__ = _make_rop(operator.truediv)

    __floordiv__ = __truediv__
    __rfloordiv__ = __rtruediv__

    def __neg__(self):
        return Vector(*-self.v)

    def __getattr__(self, symbol):
        return self.v[self.sym_to_index[symbol]]

    def __iter__(self):
        return iter(self.v)

    def length(self):
        return np.linalg.norm(self.v)

    def normalize(self):
        return self / self.length()

    def rotate(self, angle, vec):
        vec = vec.normalize()

        x, y, z = vec
        c = math.cos(math.pi/180 * angle)
        cinv = 1 - c
        s = math.sin(math.pi/180 * angle)

        matrix = np.array([[x**2*cinv + c,  x*y*cinv - z*s, x*z*cinv + y*s],
                           [y*x*cinv + z*s, y**2*cinv + c,  y*z*cinv - x*s],
                           [x*z*cinv - y*s, y*z*cinv + x*s, z**2*cinv + c ]])

        result = np.dot(matrix, self.v)
        return Vector(*result)

    def __repr__(self):
        return 'Vector(' + ', '.join(map(repr, self.v)) + ')'

    def angle_to(self, other):
        c = self.normalize() * other.normalize()

        # Occasionally, identical or almost identical
        # normalized vectors can have dot products slightly outside
        # the expected range due to floating point funness.
        # Correct that by rounding to 1 or -1 when outside
        # the range.
        if c > 1:
            c = 1
        elif c < -1:
            c = -1

        return (180/math.pi) * math.acos(c)

    def right_hand_rule_angle_to(self, other, axis):
        dot = self * other
        det = axis.normalize() * (self ** other)

        return (180/math.pi) * math.atan2(det, dot)
