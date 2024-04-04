class CubicPolynome:
    def __init__(self, a: float, b: float, c: float, d: float, x0: float) -> None:
        self.a: float = a
        self.b: float = b
        self.c: float = c
        self.d: float = d
        self.x0: float = x0

    def get(self, x) -> float:
        xoff = x - self.x0
        return self.a + self.b * xoff + self.c * (xoff**2) + self.d * (xoff**3)
