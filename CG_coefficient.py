from scipy import linalg
import numpy as np
from math import sqrt

def CG_coefficients(j_1, m_1, j_2, m_2, J, M):
    assert j_1 >= j_2, 'j_1 should not be smaller than j_2'
    if J > j_1 + j_2 or J < j_1 - j_2:
        return 0
    elif J == j_1 + j_2:
        if m_1 + m_2 != M or m_1 > j_1 or m_1 < -j_1 or m_2 > j_2 or m_2 < -j_2:
            return 0
        elif M == J:
            return 1
        else:
            return sqrt((j_1*(j_1+1)-m_1*(m_1+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1+1, j_2, m_2, J, M+1) + sqrt((j_2*(j_2+1)-m_2*(m_2+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1, j_2, m_2+1, J, M+1)
    else:
        if m_1 + m_2 != M or m_1 > j_1 or m_1 < -j_1 or m_2 > j_2 or m_2 < -j_2:
            return 0
        elif M == J:
            summary = sum([CG_coefficients(j_1, j_1, j_2, M-j_1, j, M) ** 2 for j in range(J, j_1 + j_2 + 1)])
            if m_1 == j_1:
                return (1 - summary) ** 0.5
            else:
                def solution(x):
                    A = np.array([[CG_coefficients(j_1, j_1-i, j_2, M-j_1+i, j, M) for i in range(1, j_2+1)] for j in range(j_1+j_2, J, -1)])
                    b = np.array([-CG_coefficients(j_1, j_1, j_2, M-j_1, j, M) * CG_coefficients(j_1, j_1, j_2, M-j_1, J, M) for j in range(j_1+j_2, J, -1)])
                    sol = linalg.solve(A, b)
                    return sol[x]
                return solution(j_1 - m_1 - 1)
        else:
            return sqrt((j_1*(j_1+1)-m_1*(m_1+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1+1, j_2, m_2, J, M+1) + sqrt((j_2*(j_2+1)-m_2*(m_2+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1, j_2, m_2+1, J, M+1)