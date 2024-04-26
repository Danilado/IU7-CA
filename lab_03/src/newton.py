from .matrix3d import Matrix3d
from .matrix4d import Matrix4d
from .mat_utils import get_closest_in_mat


def newton_divided_difference(mat):
    diff_mat = [[row[1] for row in mat]]
    for i in range(1, len(mat)):
        tmp = []
        for j in range(len(mat) - i):
            tmp.append(
                (diff_mat[i - 1][j + 1] - diff_mat[i - 1][j]) / (mat[i + j][0] - mat[j][0]))
        diff_mat.append(tmp)
    return diff_mat


def newton(mat, n, x):
    new_mat = get_closest_in_mat(mat, n + 1, x)
    diff_table = newton_divided_difference(new_mat)
    res = diff_table[0][0]
    x_mult = 1
    for i in range(1, n + 1):
        x_mult *= (x - new_mat[i - 1][0])
        res += x_mult * diff_table[i][0]
    return res


def newtonSecondDeriativeThird(mat, x):
    new_mat = get_closest_in_mat(mat, 4, x)
    diff_table = newton_divided_difference(new_mat)
    return 2 * diff_table[2][0] + diff_table[3][0] * (6 * x - 2 * (new_mat[0][0] + new_mat[1][0] + new_mat[2][0]))


def newton3d(mat: Matrix3d, nx, ny, x, y):
    neary = mat.get_closest_y_index(y, ny + 1)
    neary, nearyindexes = zip(*neary)
    zy = []
    for yi in nearyindexes:
        zy.append(newton(mat.get_zx_from_y(yi), nx, x))
    return newton(list(zip(neary, zy)), ny, y)


def newton4d(mat: Matrix4d, nx, ny, nz, x, y, z):
    # nearz = [i[0] for i in get_closest_in_mat([[z] for z in mat.z], nz + 1, z)]
    nearz = mat.get_closest_z_index(z, nz + 1)
    nearz, nearzindexes = zip(*nearz)
    fz = []
    for zi in nearzindexes:
        fz.append(newton3d(mat.get_z_slice(zi), nx, ny, x, y))
    return newton(list(zip(nearz, fz)), nz, z)
