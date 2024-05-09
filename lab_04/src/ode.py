from .gauss import Gauss


def ODE(x, n, funcs, coeffFuncs):
    if (len(funcs) != n + 1 or len(coeffFuncs) != n + 1):
        raise ValueError("Number of functions must equal n + 1")

    f0 = coeffFuncs[0]
    cCoeffs = [[sum(list(f1(x) * f2(x)))
                for f2 in coeffFuncs[1::]] for f1 in coeffFuncs[1::]]
    fCoeffs = [-sum(list(f0(x) * f(x))) for f in coeffFuncs[1::]]


    for c, f in zip(cCoeffs, fCoeffs):
        for coeff in c:
            print(f"{coeff:10.6f}", end="  ")
        print(f"|  {f:10.6f}")


    poliCoeffs = Gauss(cCoeffs, fCoeffs)

    print(f"solution:")
    for coeff in poliCoeffs:
        print(f"{coeff:10.6f}", end="  ")
    print()
    print()
    

    def func(xp):
        return funcs[0](xp) + sum([poliCoeffs[i] * funcs[i + 1](xp) for i in range(n)])

    return func, poliCoeffs
