from . import newton
from . import spline
from .matrix3d import Matrix3d
from .matrix4d import Matrix4d

from typing import Literal


def customInterpolation3d(
    mat: Matrix3d,
    nx: int,
    ny: int,
    x: float,
    y: float,
    x_method: Literal["newton", "spline"],
    y_method: Literal["newton", "spline"],
) -> float:
    if x_method not in ["newton", "spline"] or y_method not in ["newton", "spline"]:
        raise ValueError("Methods can only be 'newton' or 'spline'")
    zy = []
    # Если ньютон по y, то нам нужно аппроксимировать по x только в ближаших к y точках
    if (y_method == "newton"):
        neary = mat.get_closest_y_index(y, ny + 1)
        neary, nearyindexes = zip(*neary)
    else:
        neary, nearyindexes = mat.y, list(range(len(mat.y)))

    # Аппроксимуруем по X во всех точках neary со значением x
    if (x_method == "newton"):
        for yi in nearyindexes:
            zy.append(newton.newton(mat.get_zx_from_y(yi), nx, x))
    else:
        for yi in nearyindexes:
            sp = spline.Spline(mat.get_zx_from_y(yi), 0, 0)
            zy.append(sp.approximate(x))
    # Аппроксимируем по y
    if (y_method == "newton"):
        return (newton.newton(list(zip(neary, zy)), ny, y))
    else:
        sp = spline.Spline(list(zip(neary, zy)), 0, 0)
        return sp.approximate(y)


def customInterpolation4d(
        mat: Matrix4d,
        nx: int,
        ny: int,
        nz: int,
        x: float,
        y: float,
        z: float,
        x_method: Literal["newton", "spline"],
        y_method: Literal["newton", "spline"],
        z_method: Literal["newton", "spline"]
) -> float:
    if x_method not in ["newton", "spline"] or y_method not in ["newton", "spline"] or z_method not in ["newton", "spline"]:
        raise ValueError("Methods can only be 'newton' or 'spline'")
    fz = []
    if (z_method == "newton"):
        nearz = mat.get_closest_z_index(z, nz + 1)
        nearz, nearzindexes = zip(*nearz)
        for zi in nearzindexes:
            fz.append(customInterpolation3d(mat.get_z_slice(
                zi), nx, ny, x, y, x_method, y_method))
        return newton.newton(list(zip(nearz, fz)), nz, z)
    else:
        for zi in range(len(mat.z)):
            fz.append(customInterpolation3d(
                mat.f[zi], nx, ny, x, y, x_method, y_method))
        sp = spline.Spline(list(zip(mat.z, fz)), 0, 0)
        return sp.approximate(z)
