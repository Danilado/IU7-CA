from typing import Callable
from numpy import linspace


def SimpsonIntegral(start: float | int, end: float | int, deg: int, f: Callable):
    x = list(linspace(start, end, deg))

    def simfunc(a, b):
        return ((b - a) / 6) * (f(a) + 4 * f((a + b) / 2) + f(b))

    return sum([simfunc(x[i], x[i + 1]) for i in range(len(x) - 1)])
