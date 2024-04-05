import numpy as np
from . import points

Points = points.Points


def newton_2nd_deg(pts: Points, x: float) -> float:
    _, ptprev = pts.find_prev_point(x)
    data = pts.get_interval_slice(ptprev.x, 2 + 1)

    coeff_table: list[list[float]] = [
        [s.x for s in data],
        [s.y for s in data],
    ]

    for _ in range(2):
        coeff_table.append([])

    for i in range(2):
        for j in range(2 - i):
            coeff_table[i + 2].append(
                (coeff_table[i + 1][j] - coeff_table[i + 1][j + 1])
                / (coeff_table[0][j] - coeff_table[0][j + i + 1])
            )

    # print("Таблица коэф-тов:")
    # for i, line in enumerate(coeff_table):
    #     if (i < 2):
    #         print(" x:" if i == 0 else " y:", end='\t')
    #     else:
    #         print(f"y{i-2}:", end='\t')

    #     for item in line:
    #         print(f"{item:.3f}", end='\t')
    #     print()

    return 2 * coeff_table[2][0] + coeff_table[3][0] * (6 * x - 2 * (data[0].x + data[1].x + data[2].x))


def newton_interpolate(x: float, deg: int, pts: Points):
    _, ptprev = pts.find_prev_point(x)
    data = pts.get_interval_slice(ptprev.x, deg + 1)

    coeff_table: list[list[float]] = [
        [s.x for s in data],
        [s.y for s in data],
    ]

    for _ in range(deg):
        coeff_table.append([])

    for i in range(deg):
        for j in range(deg - i):
            coeff_table[i + 2].append(
                (coeff_table[i + 1][j] - coeff_table[i + 1][j + 1])
                / (coeff_table[0][j] - coeff_table[0][j + i + 1])
            )

    if False:
        print("Таблица коэф-тов:")
        for i, line in enumerate(coeff_table):
            if (i < 2):
                print(" x:" if i == 0 else " y:", end='\t')
            else:
                print(f"y{i-2}:", end='\t')

            for item in line:
                print(f"{item:.3f}", end='\t')
            print()

    res = coeff_table[1][0]
    xcoef: float = 1.0

    for i in range(deg):
        xcoef *= x - coeff_table[0][i]
        res += coeff_table[i + 2][0] * xcoef

    return res
