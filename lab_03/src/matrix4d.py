from .matrix3d import Matrix3d
from .mat_utils import get_closest_in_mat


class Matrix4d:
    def __init__(self, fname) -> None:
        self.z = []
        self.f = []
        if (fname != None):
            with open(fname, "r") as fin:
                self.read_matrix(fin)

    def read_matrix(self, f):
        s = f.readline().strip()
        n = 1
        while s != "":
            try:
                z = float(s)
            except ValueError:
                raise ValueError(f"Incorrect z format in row №{n}")

            mat = Matrix3d(None)
            try:
                mat.read_matrix(f)
            except ValueError as e:
                raise ValueError("Error while reading xy matrix") from e
            self.z.append(z)
            self.f.append(mat)

            s = f.readline().strip()
            n += 1
        if not self.check_axis():
            raise ValueError("Different xy axis in z-slices")
        self.sort_z()

    def check_axis(self):
        for i in range(len(self.f) - 1):
            if not self.f[i].axiseq(self.f[i - 1]):
                return False
        return True

    def sort_z(self):
        self.z, self.f = zip(*sorted(zip(self.z, self.f), key=lambda t: t[0]))

    def get_closest_z_index(self, z, n):
        if len(self.z) < n:
            raise ValueError(
                f"Size out of range: Z-size({len(self.z)}), n({n})")
        t = list(zip(self.z, range(len(self.z))))
        return get_closest_in_mat(t, n, z)

    def get_z_slice(self, zindex):
        if len(self.z) <= zindex:
            raise ValueError(f"Z-slice index out of range: {zindex}")
        return self.f[zindex]

    def get_max_x(self):
        return self.f[0].x[-1]

    def get_min_x(self):
        return self.f[0].x[0]

    def get_max_y(self):
        return self.f[0].y[-1]

    def get_min_y(self):
        return self.f[0].y[0]

    def generate_matrix(self, xs, xe, xsteps, ys, ye, ysteps, zs, ze, zsteps, func):
        x = [(xe - xs) / xsteps * i + xs for i in range(xsteps + 1)]
        y = [(ye - ys) / ysteps * i + ys for i in range(ysteps + 1)]
        z = [(ze - zs) / zsteps * i + zs for i in range(zsteps + 1)]
        for zv in z:
            m = Matrix3d(None)
            m.x = x
            m.y = y
            m.z = [[0] * len(x) for _ in range(len(y))]
            for i in range(len(x)):
                for j in range(len(y)):

                    # if abs(x[i] + y[j]) < 10e-8:
                    #     f = 1e8
                    # else:
                    #     f = func(x[i], y[j], zv)
                    f = func(x[i], y[j], zv)
                    m.z[j][i] = f
            self.f.append(m)
        self.z = z

    def __str__(self) -> str:
        res = ""

        for i in range(len(self.z)):
            res += f"layer: {self.z[i]}\n"
            res += str(self.f[i])
            res += "\n"

        return res
