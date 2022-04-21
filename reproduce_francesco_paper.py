"""
This is codes to trace the analysis represented by Francesco et al., 2019 paper
for finding proper method to measure the similarity bwt sample and COSMIC mutational signature.

A method for calculating linear combination of two COSMIC signature cannot be reproduced... In paper,

"... To do so, cosine similarities between the extracted signatures and each COSMIC signature,
 or a linear combination of two COSMIC signatures 
 (using non-negative least squares R package NNLS), were computed. ..."

I tried to solve linear combination of two COSMIC signature with NNLS package,
but constraints where sum of coefficients be 1.0 and both of coefficients be positive
make impossible to solve problem... NNLS package normally makes 'intercept' term.
"""

from typing import Tuple
import pandas as pd
import numpy as np
from numpy.linalg import norm
from itertools import combinations
import cvxpy as cp


class SolveLinearCombination:
    """This is class for solving below problem.

        argmin_x || A * x - y ||_2

        ( where A = (96, 2), y = (96,),
            sum(x) = 1 and all elements of x >= 0 )

    'A' matrix represents one pair of COSMIC reference mutational signatures.
    'y' vector represents a sample's mutational signature.
    'x' vector represents coefficients for pair of COSMIC signatures.
    """

    def __init__(self,
        A : pd.DataFrame,
        y : pd.DataFrame
        ) -> None:

        self._A: np.array = np.array(A)
        self._y: np.array = np.array(y)
        
        return
    
    def solve_with_sum_square(self,
        A : np.array,
        y : np.array
        ) -> Tuple[float, np.array]:

        x = cp.Variable(A.T.shape[0])
        cost = cp.sum_squares(A @ x - y)
        prob = cp.Problem(cp.Minimize(cost), [x >= 0, sum(x) == 1.0])
        prob.solve()

        return (prob.value, x.value)

    def calculate_cosine_similarity(self,
        x : np.array,
        y : np.array
        ) -> float:

        return np.dot(x, y) / (norm(x) * norm(y))