from typing import Callable


def HalfDivision(
    func: Callable[..., float | int],
    target: float | int,
    start: float | int,
    end: float | int,
    iter_limit: int = 10,
    eps: float = 1e-6,
):
    for iter in range(1, iter_limit + 1):
        center = (start + end) / 2

        # Если нашли искомую точку
        if abs(func(center) - target) < eps or iter == iter_limit:
            return center, iter

        # Умножаем на func(start) - tarrget, чтобы наклон функции
        # не влиял на результат
        if (func(start) - target) * (func(center) - target) < 0:
            end = center
        else:
            start = center
