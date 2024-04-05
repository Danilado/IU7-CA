from . import points
from . import cubic_poly

Points = points.Points
CubicPolynome = cubic_poly.CubicPolynome


class CubicSpline:
    def __init__(self, c_0: float = 0, c_n: float = 0) -> None:
        self.c_0 = c_0
        self.c_n = c_n
        self.a = []
        self.b = []
        self.c = []
        self.d = []

    def calc_a(self, pts: Points):
        res = []
        for pt in pts.data:
            res.append(pt.y)

        self.a = res[:-1]

    def calc_c(self, pts: Points):
        data = pts.data
        n = len(data)
        c = [0 for _ in range(n - 1)]
        c[0] = self.c_0 / 2
        # c[-1] = self.c_n / 2

        ksiarr = [0, 0]
        thetaarr = [0, self.c_0 / 2]

        for i in range(2, n):
            h1 = data[i].x - data[i - 1].x
            h2 = data[i - 1].x - data[i - 2].x

            ksiarr.append(ksi(ksiarr[-1], h1, h2))
            dy = data[i].y - data[i - 1].y
            dy1 = data[i - 1].y - data[i - 2].y
            thetaarr.append(theta(dy, dy1, h1, h2, thetaarr[-1], ksiarr[-2]))

        # print("ksi:", *ksiarr, sep="\t")
        # print("theta:", *thetaarr, sep="\t")

        c[-1] = thetaarr[-1] + (self.c_n / 2) * ksiarr[-1]

        for i in range(n - 2, 0, -1):
            c[i - 1] = ksiarr[i] * c[i] + thetaarr[i]

        self.c = c

    def calc_b(self, pts: Points):
        if (len(self.c) == 0):
            self.calc_c()

        data = pts.data
        n = len(data)
        b = []

        for i in range(n - 2):
            h = data[i + 1].x - data[i].x

            b_cur = (data[i + 1].y - data[i].y) / \
                h - (h * (self.c[i + 1] + 2 * self.c[i])) / 3

            b.append(b_cur)

        h = data[n - 1].x - data[n - 2].x
        b.append(
            (data[-1].y - data[- 2].y) / h - (h * (self.c_n/2 + 2 * self.c[-1])) / 3)

        self.b = b

    def calc_d(self, pts: Points):
        if (len(self.c) == 0):
            self.calc_c()

        data = pts.data
        d = []

        n = len(data)

        for i in range(n - 2):
            h = data[i + 1].x - data[i].x

            d_cur = (self.c[i + 1] - self.c[i]) / (3 * h)

            d.append(d_cur)

        h = data[-1].x - data[-2].x
        d.append(((self.c_n / 2) - self.c[-1]) / (3 * h))

        self.d = d

    def calc_coeffs(self, pts: Points):
        self.calc_a(pts)
        self.calc_c(pts)
        self.calc_b(pts)
        self.calc_d(pts)

        # print("a:", *self.a, sep="\t")
        # print("b:", *self.b, sep="\t")
        # print("c:", *self.c, sep="\t")
        # print("d:", *self.d, sep="\t")

        # print(len(self.a))
        # print(len(self.b))
        # print(len(self.c))
        # print(len(self.d))

    def get_cubic_poly(self, pts: Points, x: float) -> CubicPolynome:
        if len(self.c) == 0:
            self.calc_coeffs(pts)

        pos = 0
        while (pos < len(pts.data) - 2 and pts.data[pos+1].x < x):
            pos += 1

        # print(i)
        return CubicPolynome(
            self.a[pos],
            self.b[pos],
            self.c[pos],
            self.d[pos],
            pts.data[pos].x,)


def ksi(ksi1, h1, h2):
    return - h1 / (2 * (h1 + h2) + h2 * ksi1)


def theta(dy1, dy2, h1, h2, teta1, ksi1):
    return (3 * (dy1/h1 - dy2/h2) - h2 * teta1) / (h2 * ksi1 + 2 * (h1 + h2))
