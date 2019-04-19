import olll
import olll_numpy
import pytest


def test_readme_case():

    reduced_basis = olll.reduction([
      [1, 1, 1],
      [-1, 0, 2],
      [3, 5, 6],
    ], 0.75)

    assert reduced_basis == [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]

def test_readme_case_numpy():

    reduced_basis = olll_numpy.reduction([
      [1, 1, 1],
      [-1, 0, 2],
      [3, 5, 6],
    ], 0.75)

    assert reduced_basis == [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
