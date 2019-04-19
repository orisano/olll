from typing import Sequence, List

import numpy as np
from fractions import Fraction


class Vector(object):
    def __init__(self, x):
        fractions = list(map(Fraction, x))
        self._numerators = np.array(list(map(lambda x: x.numerator, fractions)))
        self._denominators = np.array(list(map(lambda x: x.denominator, fractions)))

    def sdot(self) -> Fraction:
        return self.dot(self)

    def dot(self, rhs: "Vector") -> Fraction:
        assert len(self) == len(rhs)
        return Fraction(np.dot(self._numerators, rhs._numerators), np.dot(self._denominators, rhs._denominators))

    def proj_coff(self, rhs: "Vector") -> Fraction:
        """
        >>> Vector([1, 1, 1]).proj_coff([-1, 0, 2])
        Fraction(1, 3)
        """
        assert len(self) == len(rhs)
        return Fraction(self.dot(rhs), self.sdot())

    def proj(self, rhs: "Vector") -> "Vector":
        """
        >>> Vector([1, 1, 1]).proj([-1, 0, 2])
        [1/3, 1/3, 1/3]
        """
        assert len(self) == len(rhs)
        return self.proj_coff(rhs) * self

    def __sub__(self, rhs: "Vector") -> "Vector":
        """
        >>> Vector([1, 2, 3]) - [6, 5, 4]
        [-5, -3, -1]
        """
        assert len(self) == len(rhs)

        return Vector._create_from_numerators_and_dominators(
            self._numerators * rhs._denominators - rhs._numerators * self._denominators,
            self._denominators * rhs._denominators
        )

    def __mul__(self, rhs: Fraction) -> "Vector":
        """
        >>> Vector(["3/2", "4/5", "1/4"]) * 2
        [3, 8/5, 1/2]
        """
        return Vector._create_from_numerators_and_dominators(
            self._numerators * rhs.numerator,
            self._denominators * rhs.denominator,
        )

    def __rmul__(self, lhs: Fraction) -> "Vector":
        """
        >>> 2 * Vector(["3/2", "4/5", "1/4"])
        [3, 8/5, 1/2]
        """
        return Vector._create_from_numerators_and_dominators(
            lhs.numerator * self._numerators,
            lhs.denominator * self._denominators,
        )

    def __len__(self):
        return self._numerators.shape[0]

    @classmethod
    def _create_from_numerators_and_dominators(cls, numerators, denominators):
        v = Vector([])
        v._denominators = denominators
        v._numerators = numerators

        return v

    def any(self):
        return self._numerators.any()

    def get_integers(self):
        assert abs(self._numerators % self._denominators).sum() == 0
        return list(self._numerators // self._denominators)


def gramschmidt(v: Sequence[Vector]) -> Sequence[Vector]:
    """
    >>> gramschmidt([[3, 1], [2, 2]])
    [[3, 1], [-2/5, 6/5]]
    >>> gramschmidt([[4, 1, 2], [4, 7, 2], [3, 1, 7]])
    [[4, 1, 2], [-8/7, 40/7, -4/7], [-11/5, 0, 22/5]]
    """
    u: List[Vector] = []
    for vi in v:
        ui = vi
        for uj in u:
            ui = ui - uj.proj(vi)

        if ui.any():
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

    return [b.get_integers() for b in basis]
