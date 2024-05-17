from typing import Callable, Optional
from .least_square import LeastSquare
from math import exp, log, sqrt


def integrate_2d(
        x_values: list[float],
        y_values: list[float],
        z_values: list[list[float]],
        limits_checker: Callable[[float, float], bool],
        x_integral_f: Callable,
        y_integral_f: Callable,
        deg_x: int,
        deg_y: int
):
    filtered_points = [[] for _ in range(len(z_values))]

    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            if (limits_checker(x, y)):
                filtered_points[j].append(log(z_values[j][i]))

    filtered_points = [*filter(lambda line: len(line) > 1, filtered_points)]

    print("matrix after filtering and aligning (with log):")
    print(r"y\x", end="\t")
    print(*[f"{x:+6.5f}" for x in x_values], sep="  ")
    for y, line in zip(y_values, filtered_points):
        print(f"{y:+4.2f}", end="\t")
        print(*[f"{z:+6.5f}" for z in line], sep="  ")
    print()

    x_integrals: list[float] = []
    for y, line in zip(y_values, filtered_points):
        if (len(line) <= 1):
            x_integrals.append(0)
            continue

        x_line = x_values[0:len(line)]

        _f1 = LeastSquare(x_line, line, [], 4)

        def f1(x: float) -> float:
            return exp(_f1(x))

        x_start = x_line[0]
        x_end = x_line[-1]

        x_integrals.append(x_integral_f(x_start, x_end, deg_x, f1))

    y_line = y_values[0:len(x_integrals)]

    print("y values        : ", end="")
    print(*[f"{y:+4.2f}" for y in y_line], sep="\t\t")
    print("integral values : ", end="")
    print(*[f"{integral:+10.8f}" for integral in x_integrals], sep="\t")
    print()

    f2 = LeastSquare(y_line,
                     x_integrals, [], 4)

    y_min = y_line[0]
    y_max = y_line[-1]

    return y_integral_f(y_min, y_max, deg_y, f2)


def integrate_2d_concrete(
        x_values: list[float],
        y_values: list[float],
        concrete_function: Callable,
        limits_checker: Callable[[float, float], bool],
        x_limit_getter: Callable,
        x_integral_f: Callable,
        y_integral_f: Callable,
        deg_x: int,
        deg_y: int,
):
    x_integrals: list[float] = []
    y_line = []
    for y in y_values:
        x_start = x_limit_getter(y)
        if (x_start == sqrt(y)):
            x_end = sqrt(y)
            x_start = - sqrt(y)
        else:
            x_end = x_values[-1]

        if not limits_checker(x_start, y):
            continue

        def f1(x: float) -> float:
            return concrete_function(x, y)

        if abs(x_end - x_start) < 1e-6:
            continue

        y_line.append(y)
        x_integrals.append(x_integral_f(x_start, x_end, deg_x, f1))

    print("y values        : ", end="")
    print(*[f"{y:+4.2f}" for y in y_line], sep="\t\t")
    print("integral values : ", end="")
    print(*[f"{integral:+10.8f}" for integral in x_integrals], sep="\t")
    print()

    f2 = LeastSquare(y_line,
                     x_integrals, [], len(x_integrals))

    y_min = y_line[0]
    y_max = y_line[-1]

    return y_integral_f(y_min, y_max, deg_y, f2)
