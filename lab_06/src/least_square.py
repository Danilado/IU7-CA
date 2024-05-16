def __Gauss(A, B):
    """
        A | B
        A -> triangle
        A -> diagonal

        res[i] = A[i][i] / B[i] 
    """

    n = len(A)
    if n != len(B):
        raise ValueError("Mat and Rights must have the same length")

    # приводим к треугольному виду
    for i in range(n):
        for j in range(i + 1, n):
            coeff = -(A[j][i] / A[i][i])
            for k in range(i, n):
                A[j][k] += coeff * A[i][k]
            B[j] += coeff * B[i]

    # Приводим к диагональному виду
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            coeff = -(A[j][i] / A[i][i])
            A[j][i] += coeff * A[i][i]
            B[j] += coeff * B[i]

    #  Находим решения слау
    res = [B[i] / A[i][i] for i in range(n)]
    return res


def __LeastSquaresCoeffs(x, y, p, n):
    dataLen = len(x)
    if dataLen != len(y) or dataLen != len(p):
        raise ValueError("x, y and p must have the same length")

    aCoeffsx = [sum([p[i] * x[i] ** j for i in range(dataLen)])
                for j in range(2 * n + 1)]
    aCoeffsy = [sum([p[i] * y[i] * x[i] ** j for i in range(dataLen)])
                for j in range(n + 1)]
    CoeffMat = [[aCoeffsx[i + j] for i in range(n + 1)] for j in range(n + 1)]

    return __Gauss(CoeffMat, aCoeffsy)


def LeastSquare(x, y, weights, deg):
    if (weights == []):
        weights = [1 for _ in range(len(x))]
    a = __LeastSquaresCoeffs(x, y, weights, deg)

    def func(xp):
        return sum([a[i] * xp ** i for i in range(deg + 1)])

    return func
