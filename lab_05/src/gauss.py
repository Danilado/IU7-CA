def Gauss(A, B):
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
