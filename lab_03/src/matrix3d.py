from .mat_utils import get_closest_in_mat

EPS = 1e-3


class Matrix3d:
    def __init__(self, fname) -> None:
        self.x = []
        self.y = []
        self.z = []
        if (fname != None):
            with open(fname, "r") as fin:
                self.read_matrix(fin)

    def read_matrix(self, fin):
        x = fin.readline().strip().split()
        if len(x) <= 1:
            raise ValueError("Empty x-axis")

        x = x[1:]
        try:
            self.x = list(map(float, x))
        except ValueError:
            raise ValueError(f"Incorrect matrix x-axis format:{' '.join(x)}")
        n = 1
        self.y = []
        self.z = []
        s = fin.readline().strip()
        while s != "":
            s = s.split()
            if len(s) - 1 != len(self.x):
                raise ValueError(
                    f"Incorrect row №{n} size: x-axis:{len(self.x)}, row:{len(s) - 1}")

            try:
                s = list(map(float, s))
            except ValueError:
                raise ValueError(
                    f"Incorrect matrix row №{n} format:{' '.join(x)}")

            self.y.append(s[0])
            self.z.append(s[1:])
            s = fin.readline().strip()
            n += 1

        self.sort_x()
        self.sort_y()

    def sort_x(self):
        self.transpose_z()
        self.x, self.z = zip(*sorted(zip(self.x, self.z), key=lambda t: t[0]))
        self.transpose_z()

    def sort_y(self):
        self.y, self.z = zip(*sorted(zip(self.y, self.z), key=lambda t: t[0]))

    def transpose_z(self):
        tz = [[0] * len(self.x) for _ in range(len(self.y))]
        for i in range(len(self.y)):
            for j in range(len(self.x)):
                tz[i][j] = self.z[j][i]
        self.z = tz

    def get_zx_from_y(self, yindex):
        mat = []
        for i in range(len(self.x)):
            mat.append([self.x[i], self.z[yindex][i]])
        return mat

    def get_closest_y_index(self, y, n):
        if len(self.y) < n:
            raise ValueError(
                f"Size out of range: Y-size({len(self.y)}), n({n})")
        t = list(zip(self.y, range(len(self.y))))
        return get_closest_in_mat(t, n, y)

    def axiseq(self, __value: object) -> bool:
        if len(self.x) != len(__value.x):
            return False
        if len(self.y) != len(__value.y):
            return False

        for i in range(len(self.x)):
            if abs(self.x[i] - __value.x[i]) > EPS:
                return False

        for i in range(len(self.y)):
            if abs(self.y[i] - __value.y[i]) > EPS:
                return False

        return True

    def __str__(self) -> str:
        res = ""
        res += (f"\tx: {self.x}\n")
        res += (f"\ty: {self.y}\n")
        res += (f"\tz:\n")
        for z in self.z:
            res += (f"\t\t")
            for zz in z:
                res += (f"{zz:0.3f} ")
            res += (f"\n")
        return res
