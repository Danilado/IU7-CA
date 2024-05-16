from typing import Callable, Tuple
from math import cos, pi
from functools import cache


EPS = 1e-16


def __LegendreEvaluate(x: float, deg: int) -> Tuple[float, float]:
    vbuf1 = x
    vbuf2 = 1

    fract = 1 / (x**2 - 1)

    for i in range(2, deg+1):
        v = ((2*i - 1) * x * vbuf1 - (i-1) * vbuf2) / i
        d = i * fract * (x * v - vbuf1)

        vbuf2 = vbuf1
        vbuf1 = v

    return v, d


def print_legendre(deg: int):
    roots, weights = LegendrePolynome(deg)
    print(f"for deg: {deg}")
    print("\nLegendre roots:")
    for r in roots:
        print(f"{r:+10.8f}", end="  ")
    print("\nLegendre weights:")
    for w in weights:
        print(f"{w:+10.8f}", end="  ")
    print()


@cache
def LegendrePolynome(deg: int) -> Tuple[list[float], list[float]]:
    if deg == 0:
        return [0, 0], [0, 2]

    if deg == 1:
        return [0, +0.57735, - 0.57735], [0, 1, 1]

    roots: list[float] = [0 for _ in range(deg+1)]
    weights: list[float] = [0 for _ in range(deg+1)]
    for i in range(deg+1):
        x = cos(pi * (i - 0.25) / (deg + 0.5))
        v, d = __LegendreEvaluate(x, deg)
        dr = v / d

        while (abs(dr) > EPS):
            x = x - dr
            v, d = __LegendreEvaluate(x, deg)
            dr = v / d

        roots[i] = x
        weights[i] = 2 / ((1 - x**2) * d**2)

    roots = roots[1::]
    roots = roots[::-1]
    weights = weights[1::]
    weights = weights[::-1]

    return roots, weights


def GaussIntegral(start: float | int, end: float | int, deg: int, f: Callable):
    p = (end - start)/2
    q = (start + end)/2

    roots, weights = LegendrePolynome(deg)

    if False:
        print("\nLegendre roots:")
        for r in roots[1::]:
            print(f"{r:+8.6f}", end="  ")
        print("\nLegendre weights:")
        for w in weights[1::]:
            print(f"{w:+8.6f}", end="  ")
        print()

    s: float = 0.0
    for i in range(deg):
        s += weights[i] * f(p * roots[i] + q)

    return p * s
