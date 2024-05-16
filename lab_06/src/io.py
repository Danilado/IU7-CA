from typing import Tuple


def read_3d(filename: str) -> Tuple[list[float], list[float], list[list[float]]]:
    with open(filename, "r") as f:
        x: list[float]
        y: list[float] = []
        z: list[list[float]]
        x = list(map(float, f.readline().strip().split()[1::]))
        other = f.readlines()

        z = [[] for _ in range(len(other))]
        for i, line in enumerate(other):
            l = list(map(float, line.split()))
            y.append(l.pop(0))
            z[i] = [*l]

        return x, y, z
