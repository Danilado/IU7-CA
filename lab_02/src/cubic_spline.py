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
        c = [0 for _ in range(n)]
        c[0] = self.c_0 / 2
        c[-1] = self.c_n / 2

        ksi = [0, 0]
        theta = [0, 0]

        for i in range(2, n):
            h1 = data[i].x - data[i - 1].x
            h2 = data[i - 1].x - data[i - 2].x

            phi = 3 * ((data[i].y - data[i - 1].y) /
                       h1 - (data[i - 1].y - data[i - 2].y) / h2)

            ksi_cur = - h1 / (h2 * ksi[i - 1] + 2 * (h2 + h1))
            theta_cur = (phi - h1 * theta[i - 1]) / \
                (h1 * ksi[i - 1] + 2 * (h2 + h1))

            ksi.append(ksi_cur)
            theta.append(theta_cur)

        c[-2] = theta[len(theta) - 1]

        for i in range(n - 2, 0, -1):
            c[i - 1] = ksi[i] * c[i] + theta[i]

        self.c = c

    def calc_b(self, pts: Points):
        if (len(self.c) == 0):
            self.calc_c()

        data = pts.data
        n = len(data)
        b = []

        for i in range(1, n - 1):
            h = data[i].x - data[i - 1].x

            b_cur = (data[i].y - data[i - 1].y) / \
                h - (h * (self.c[i] + 2 * self.c[i - 1])) / 3

            b.append(b_cur)

        h = data[n - 1].x - data[n - 2].x
        b.append(
            (data[n - 1].y - data[n - 2].y) / h - (h * 2 * self.c[i]) / 3)

        self.b = b

    def calc_d(self, pts: Points):
        if (len(self.c) == 0):
            self.calc_c()

        data = pts.data
        d = []

        n = len(data)

        for i in range(1, n - 1):
            h = data[i].x - data[i - 1].x

            d_cur = (self.c[i] - self.c[i - 1]) / (3 * h)

            d.append(d_cur)

        h = data[n - 1].x - data[n - 2].x
        d.append((- self.c[i]) / (3 * h))

        self.d = d

    def calc_coeffs(self, pts: Points):
        self.calc_a(pts)
        self.calc_c(pts)
        self.calc_b(pts)
        self.calc_d(pts)

        print("a:", *self.a, sep="\t")
        print("b:", *self.b, sep="\t")
        print("c:", *self.c, sep="\t")
        print("d:", *self.d, sep="\t")

    def get_cubic_poly(self, pts: Points, x: float) -> CubicPolynome:
        if len(self.c) == 0:
            self.calc_coeffs(pts)

        pos = 1
        while (pos < len(pts.data) and pts.data[pos].x < x):
            pos += 1
        pos -= 1

        # print(i)
        return CubicPolynome(self.a[pos], self.b[pos], self.c[pos], self.d[pos], pts.data[pos].x)
