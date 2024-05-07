from .gauss import Gauss


def LeastSquaresCoeffs(x, y, p, n):
    dataLen = len(x)
    if dataLen != len(y) or dataLen != len(p):
        raise ValueError("x, y and p must have the same length")

    aCoeffsx = [sum([p[i] * x[i] ** j for i in range(dataLen)])
                for j in range(2 * n + 1)]
    aCoeffsy = [sum([p[i] * y[i] * x[i] ** j for i in range(dataLen)])
                for j in range(n + 1)]
    CoeffMat = [[aCoeffsx[i + j] for i in range(n + 1)] for j in range(n + 1)]

    return Gauss(CoeffMat, aCoeffsy)


def LeastSquare(x, y, p, n):
    a = LeastSquaresCoeffs(x, y, p, n)

    def func(xp):
        return sum([a[i] * xp ** i for i in range(n + 1)])

    return func


def LeastSquaresCoeffs2D(x, y, z, p, n):
    dataLen = len(x)
    if dataLen != len(y) or dataLen != len(p) or dataLen != len(p):
        raise ValueError("x, y, z and p must have the same length")

    def aSum(xdeg, ydeg):
        return sum([xn ** xdeg * yn ** ydeg * pn for xn, yn, zn, pn in zip(x, y, z, p)])

    def zSum(xdeg, ydeg):
        return sum([xn ** xdeg * yn ** ydeg * zn * pn for xn, yn, zn, pn in zip(x, y, z, p)])

    aCoeffs = [[aSum(i + k, j + t) for k in range(n + 1) for t in range(n + 1 - k)]
               for i in range(n + 1) for j in range(n + 1 - i)]
    zCoeffs = [zSum(i, j) for i in range(n + 1) for j in range(n + 1 - i)]

    return Gauss(aCoeffs, zCoeffs)


def LeastSquare2D(x, y, z, p, n):
    a = LeastSquaresCoeffs2D(x, y, z, p, n)

    def func(xp, yp):
        res = 0
        ind = 0
        for i in range(n + 1):
            for j in range(n + 1 - i):
                res += a[ind] * xp ** i * yp ** j
                ind += 1
        return res
    return func


def ScalarMult(x, f1, f2):
    return sum([f1(xn) * f2(xn) for xn in x])


def yScalarMult(x, y, f):
    return sum([f(x[i]) * y[i] for i in range(len(x))])
