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
            if m_1 == j_1:
                summary = sum([CG_coefficients(j_1, j_1, j_2, M-j_1, j/2, M) ** 2 for j in range(int(2*J+2), int(2*(j_1 + j_2 + 1)), 2)])
                return (1 - summary) ** 0.5
            else:
                A = np.array([[CG_coefficients(j_1, j_1-i/2, j_2, M-j_1+i/2, j/2, M) for i in range(2, int(2*(j_1+j_2-J+1)), 2)] for j in range(int(2*(j_1+j_2)), int(2*J), -2)])
                b = np.array([(-1) * CG_coefficients(j_1, j_1, j_2, M-j_1, j/2, M) * CG_coefficients(j_1, j_1, j_2, M-j_1, J, M) for j in range(int(2*(j_1+j_2)), int(2*J), -2)])
                solution = linalg.solve(A, b)
                return solution[int(j_1 - m_1 - 1)]
        else:
            return sqrt((j_1*(j_1+1)-m_1*(m_1+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1+1, j_2, m_2, J, M+1) + sqrt((j_2*(j_2+1)-m_2*(m_2+1)) / (J*(J+1)-M*(M+1))) * CG_coefficients(j_1, m_1, j_2, m_2+1, J, M+1)
