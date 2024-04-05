from . import point
Point = point.Point


class Points:
    def __init__(self, init_data: list[Point] = []) -> None:
        self.data: list[Point] = init_data

    def sort_data(self):
        self.data.sort(key=lambda p: p.x)

    def add_point(self, p: Point):
        self.data.append(p)

    def parse_table(self, t: list[list[float]]):
        for i in range(len(t[0])):
            self.data.append(Point(t[0][i], t[1][i]))
        self.sort_data()

    def find_interval(self, x: float, length: int) -> {int, int}:
        closest_i = 0, 99999999999999
        for i, p in enumerate(self.data):
            if abs(p.x - x) < closest_i[1]:
                closest_i = i, abs(p.x - x)

        i1: int
        i2: int

        if closest_i[0] >= length // 2:
            i1 = closest_i[0] - length // 2
        else:
            i1 = 0

        if len(self.data) > i1 + length:
            i2 = i1 + length
        else:
            i1 = len(self.data) - length
            i2 = i1 + length

        return i1, i2

    def get_interval_slice(self, x, length) -> list[Point]:
        i1, i2 = self.find_interval(x, length)
        return self.data[i1:i2]

    def get_dx_list(self) -> list[float]:
        res = []
        pprev = False
        for p in self.data:
            if pprev:
                res.append(p.x - pprev.x)
                pprev = p
            else:
                pprev = p

        return res

    def find_closest_point(self, x) -> {int, Point}:
        closest_i = 0, 99999999999999
        for i, p in enumerate(self.data):
            if abs(p.x - x) < closest_i[1]:
                closest_i = i, abs(p.x - x)

        return closest_i[0], self.data[closest_i[0]]

    def find_prev_point(self, x) -> {int, Point}:
        i = 0
        while i < len(self.data) - 1 and x > self.data[i].x:
            i += 1

        return i, self.data[i]

    def find_next_point(self, x) -> {int, Point}:
        i = 0
        while i < len(self.data) and x > self.data[i].x:
            i += 1

        i += 1
        if (i >= len(self.data)):
            i -= 1

        return i, self.data[i]

    def __str__(self) -> str:
        res = ""
        res += "x:\t"
        for p in self.data:
            res += f"{p.x:.3f}\t"
        res += "\n"
        res += "y:\t"
        for p in self.data:
            res += f"{p.y:.3f}\t"
        return res
