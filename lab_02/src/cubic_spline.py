from . import points
from . import cubic_poly

Points = points.Points
CubicPolynome = cubic_poly.CubicPolynome


class CubicSpline:
    def __init__(self, c_0: float = 0, c_n: float = 0) -> None:
        self.c_0 = c_0
        self.c_n = c_n
        self.ksi = []
        self.eta = []
        self.c = []

    def thomas_calc(self, pts: Points):
        h = pts.get_dx_list()
        data = pts.data

        ksi = [0]
        eta = [self.c_0 / 2]

        for i in range(1, len(h)):
            a = h[i - 1]
            b = -2 * (h[i - 1] + h[i])
            d = h[i]
            f = -3 * (
                (data[i + 1].y - data[i].y) / h[i]
                - (data[i].y - data[i - 1].y) / h[i - 1]
            )

            ksi.append(d / (b - a * ksi[i - 1]))
            eta.append((a * eta[i - 1] + f) / (b - a * ksi[i - 1]))

        n = len(data)
        c = [0 for _ in range(n)]
        c[0] = self.c_0 / 2
        c[n - 1] = self.c_n / 2

        for i in range(len(ksi) - 1, 0, -1):
            c[i] = ksi[i] * c[i + 1] + eta[i]

        self.ksi = ksi
        self.eta = eta
        self.c = c

        # print("ksi:", *self.ksi, sep="\t")
        # print("eta:", *self.eta, sep="\t")
        # print("c:", *self.c, sep="\t")

    def thomas(self, pts: Points, x: float) -> CubicPolynome:
        h = pts.get_dx_list()
        data = pts.data

        if len(self.c) == 0:
            self.thomas_calc(pts)

        i, pt = pts.find_closest_point(x)
        if i >= len(self.c) - 1:
            i = len(self.c) - 2
        return CubicPolynome(
            (self.c[i + 1] - self.c[i]) / (3.0 * h[i]),
            self.c[i + 1],
            (data[i + 1].y - data[i].y) / h[i]
            - 1.0 / 3.0 * h[i] * (self.c[i + 1] + 2.0 * self.c[i]),
            data[i].y,
            pt.x,
        )
