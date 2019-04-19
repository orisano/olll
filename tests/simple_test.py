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

    reduced_basis = olll.reduction([
        [105, 821, 404, 328],
        [881, 667, 644, 927],
        [181, 483, 87, 500],
        [893, 834, 732, 441]
    ], 0.75)

    assert [[76, -338, -317, 172], [88, -171, -229, -314], [269, 312, -142, 186], [519, -299, 470, -73]] == reduced_basis

def test_readme_case_numpy():

    reduced_basis = olll_numpy.reduction([
      [1, 1, 1],
      [-1, 0, 2],
      [3, 5, 6],
    ], 0.75)

    assert [[0, 1, 0], [1, 0, 1], [-1, 0, 2]] == reduced_basis

    reduced_basis = olll_numpy.reduction([
        [105, 821, 404, 328],
        [881, 667, 644, 927],
        [181, 483, 87, 500],
        [893, 834, 732, 441]
    ], 0.75)

    assert [[76, -338, -317, 172], [88, -171, -229, -314], [269, 312, -142, 186], [519, -299, 470, -73]] == reduced_basis
