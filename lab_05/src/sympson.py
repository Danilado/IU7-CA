from typing import Callable


def SympsonIntegral(x: list[float | int], f: Callable):

    def symfunc(a, b):
        return ((b - a) / 6) * (f(a) + 4 * f((a + b) / 2) + f(b))

    return sum([symfunc(x[i], x[i + 1]) for i in range(len(x) - 1)])
