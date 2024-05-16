def onesided(arr: list[float]) -> list[float]:
    res = []
    for i in range(len(arr)-1):
        res.append(arr[i + 1] - arr[i])
    res.append(None)

    return res


def mid(arr: list[float]) -> list[float]:
    res = [None]
    for i in range(len(arr) - 2):
        res.append((arr[i] + arr[i+1])/2)
    res.append(None)

    return res


def runge(arr: list[float]) -> list[float]:
    res = [None]

    for i in range(len(arr) - 2):
        d1 = arr[i + 1] - arr[i]
        d2 = (arr[i + 2] - arr[i])/2

        res.append(d1 + (d1 - d2) / 3)

    res.append(None)

    return res


def alignment(x_arr: list[float], y_arr: list[float]) -> list[float]:
    res = []

    for i in range(len(y_arr)-1):
        y0frac = 1/y_arr[i]
        y1frac = 1/y_arr[i+1]

        x0frac = 1/x_arr[i]
        x1frac = 1/x_arr[i+1]

        yfracdiff = y1frac - y0frac
        xfracdiff = x1frac - x0frac

        res.append((yfracdiff/xfracdiff)*(y_arr[i]**2)/(x_arr[i]**2))

    res.append(None)

    return res


def second_der(arr: list[float]) -> list[float]:
    res = [None]

    for i in range(len(arr) - 2):
        res.append(arr[i] - 2 * arr[i+1] + arr[i+2])

    res.append(None)

    return res
