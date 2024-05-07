from .gauss import Gauss


def ODE(x, n, funcs, coeffFuncs):
    if (len(funcs) != n + 1 or len(coeffFuncs) != n + 1):
        raise ValueError("Number of functions must equal n + 1")

    f0 = coeffFuncs[0]
    cCoeffs = [[sum([f1(xp) * f2(xp) for xp in x])
                for f2 in coeffFuncs[1:]] for f1 in coeffFuncs[1:]]
    fCoeffs = [-sum([f0(xp) * f(xp) for xp in x]) for f in coeffFuncs[1:]]

    print(f"for {n=}")
    print(f"cCoeffs:")
    for c in cCoeffs:
        for coeff in c:
            print(f"{coeff:10.6f}", end="  ")
        print()
    print(f"fCoeffs:")
    for coeff in fCoeffs:
        print(f"{coeff:10.6f}", end="  ")
    print()

    poliCoeffs = Gauss(cCoeffs, fCoeffs)

    print(f"solution:")
    for coeff in poliCoeffs:
        print(f"{coeff:10.6f}", end="  ")
    print()
    print()

    def func(xp):
        return funcs[0](xp) + sum([poliCoeffs[i] * funcs[i + 1](xp) for i in range(n)])

    return func, poliCoeffs
