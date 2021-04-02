from typing import Sequence, List

import numpy as np


def gramschmidt(v: np.ndarray) -> np.ndarray:
    """
    >>> gramschmidt_np(np.array([[3, 1], [2, 2]]))
    [[3, 1], [-2/5, 6/5]]
    >>> gramschmidt_np(np.array([[4, 1, 2], [4, 7, 2], [3, 1, 7]]))
    [[4, 1, 2], [-8/7, 40/7, -4/7], [-11/5, 0, 22/5]]
    """
    v = v.astype(np.float32)
    v = v.copy()
    u_sdots = np.zeros(len(v))
    u_sdots[0] = v[0].dot(v[0])
    for i in range(1, len(v)):
        vi = v[i]
        u = v[:i]
        vi = vi - (u.dot(vi) / u_sdots[:i]).sum(axis=0)
        v[i] = vi
        u_sdots[i] = vi.dot(vi)
    return v

def reduction(basis: np.ndarray, delta: float) -> Sequence[Sequence[int]]:
    """
    >>> reduction(np.array([[1, 1, 1], [-1, 0, 2], [3, 5, 6]]), 0.75)
    [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
    >>> reduction(np.array([[105, 821, 404, 328], [881, 667, 644, 927], [181, 483, 87, 500], [893, 834, 732, 441]]), 0.75)
    [[76, -338, -317, 172], [88, -171, -229, -314], [269, 312, -142, 186], [519, -299, 470, -73]]
    """
    n = len(basis)
    ortho = gramschmidt(basis)

    def mu(i: int, j: int) -> float:
        a, b = ortho[j], basis[i]
        return np.dot(a, b) / np.dot(a, a)

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            mu_kj = mu(k, j)
            if abs(mu_kj) > 0.5:
                basis[k] = basis[k] - basis[j] * round(mu_kj)
                ortho = gramschmidt(basis)

        if ortho[k].dot(ortho[k]) >= (delta - mu(k, k - 1)**2) * ortho[k - 1].dot(ortho[k - 1]):
            k += 1
        else:
            basis[k], basis[k - 1] = basis[k - 1], basis[k].copy()
            ortho = gramschmidt(basis)
            k = max(k - 1, 1)

    return basis

