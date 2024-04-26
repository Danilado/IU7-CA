def get_closest_in_mat(mat, n, x):
    if len(mat) < n:
        raise ValueError("Matrix too small")

    i = 1
    while i < len(mat) and x - mat[i][0] > 0:
        i += 1

    if i == len(mat):
        return mat[-n:]

    res = []
    j = i - 1
    while len(res) < n:
        if (i >= len(mat) and j < 0):
            raise ValueError("Matrix too small")

        if (j >= 0):
            res.append(mat[j])
            j -= 1

        if (i < len(mat) and len(res) < n):
            res.append(mat[i])
            i += 1

    return sorted(res, key=lambda x: x[0])
