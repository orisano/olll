import operator
from typing import Iterator, List, Sequence

Vector = Sequence[float]

def _dot(a: Vector, b: Vector) -> float:
    """
    >>> _dot((1.0, 2.0, 3.0), (4.0, 5.0, 6.0))
    32.0
    """
    assert len(a) == len(b)
    return sum(map(lambda x: x[0] * x[1], zip(a, b)))


def _sub(a: Vector, b: Vector) -> Vector:
    """
    >>> _sub((1.0, 2.0, 3.0), (4.0, 5.0, 6.0))
    (-3.0, -3.0, -3.0)
    """
    assert len(a) == len(b)
    return tuple(x - y for x, y in zip(a, b))


def _mul(a: Vector, b: float) -> Vector:
    """
    >>> _mul((1.0, 2.0, 3.0), 4.0)
    (4.0, 8.0, 12.0)
    """
    return tuple(x * b for x in a)


def gramschmidt(v: Sequence[Vector]) -> Iterator[Vector]:
    """
    >>> list(gramschmidt([[1, 2], [3, 4]]))
    [(1, 2), (0.8, -0.4)]
    """
    u: List[Vector] = []
    cache_dot: List[float] = []
    for vi in v:
        ui = tuple(vi)
        for uj, dot in zip(u, cache_dot):
            ui = _sub(ui, _mul(uj, _dot(uj, vi) / dot))
        yield ui
        u.append(ui)
        cache_dot.append(_dot(ui, ui))


def reduction(basis: Sequence[Sequence[int]], delta: float):
    pass
