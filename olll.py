from fractions import Fraction
from typing import List, Sequence


class Vector(list):
    def __init__(self, x):
        super().__init__(map(Fraction, x))

    def sdot(self) -> Fraction:
        return self.dot(self)

    def dot(self, rhs: "Vector") -> Fraction:
        """
        >>> Vector([1, 2, 3]).dot([4, 5, 6])
        Fraction(32, 1)
        """
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return sum(map(lambda x: x[0] * x[1], zip(self, rhs)))

    def proj_coff(self, rhs: "Vector") -> Fraction:
        """
        >>> Vector([1, 1, 1]).proj_coff([-1, 0, 2])
        Fraction(1, 3)
        """
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return self.dot(rhs) / self.sdot()

    def proj(self, rhs: "Vector") -> "Vector":
        """
        >>> Vector([1, 1, 1]).proj([-1, 0, 2])
        [1/3, 1/3, 1/3]
        """
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return self.proj_coff(rhs) * self

    def __sub__(self, rhs: "Vector") -> "Vector":
        """
        >>> Vector([1, 2, 3]) - [6, 5, 4]
        [-5, -3, -1]
        """
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return Vector(x - y for x, y in zip(self, rhs))

    def __mul__(self, rhs: Fraction) -> "Vector":
        """
        >>> Vector(["3/2", "4/5", "1/4"]) * 2
        [3, 8/5, 1/2]
        """
        return Vector(x * rhs for x in self)

    def __rmul__(self, lhs: Fraction) -> "Vector":
        """
        >>> 2 * Vector(["3/2", "4/5", "1/4"])
        [3, 8/5, 1/2]
        """
        return Vector(x * lhs for x in self)

    def __repr__(self) -> str:
        return "[{}]".format(", ".join(str(x) for x in self))


def gramschmidt(v: Sequence[Vector]) -> Sequence[Vector]:
    """
    >>> gramschmidt([[3, 1], [2, 2]])
    [[3, 1], [-2/5, 6/5]]
    >>> gramschmidt([[4, 1, 2], [4, 7, 2], [3, 1, 7]])
    [[4, 1, 2], [-8/7, 40/7, -4/7], [-11/5, 0, 22/5]]
    """
    u: List[Vector] = []
    for vi in v:
        ui = Vector(vi)
        for uj in u:
            ui = ui - uj.proj(vi)

        if any(ui):
            u.append(ui)
    return u


def reduction(basis: Sequence[Sequence[int]], delta: float) -> Sequence[Sequence[int]]:
    """
    >>> reduction([[1, 1, 1], [-1, 0, 2], [3, 5, 6]], 0.75)
    [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
    >>> reduction([[105, 821, 404, 328], [881, 667, 644, 927], [181, 483, 87, 500], [893, 834, 732, 441]], 0.75)
    [[76, -338, -317, 172], [88, -171, -229, -314], [269, 312, -142, 186], [519, -299, 470, -73]]
    """
    n = len(basis)
    basis = list(map(Vector, basis))
    ortho = gramschmidt(basis)

    def mu(i: int, j: int) -> Fraction:
        return ortho[j].proj_coff(basis[i])

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            mu_kj = mu(k, j)
            if abs(mu_kj) > 0.5:
                basis[k] = basis[k] - basis[j] * round(mu_kj)
                ortho = gramschmidt(basis)

        if ortho[k].sdot() >= (delta - mu(k, k - 1)**2) * ortho[k - 1].sdot():
            k += 1
        else:
            basis[k], basis[k - 1] = basis[k - 1], basis[k]
            ortho = gramschmidt(basis)
            k = max(k - 1, 1)

    return [list(map(int, b)) for b in basis]
