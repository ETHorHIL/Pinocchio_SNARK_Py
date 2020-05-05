import numpy as np
import code_to_r1cs
from ssbls120 import Fp, Poly, Group

text = """
def qeval(x):
    y = x**3
    return y + x + 5
"""


def testing():
    r, A, B, C = code_to_r1cs.code_to_r1cs_with_inputs(text, [3])

    L = np.array([[Fp(x) for x in A[k]] for k in range(len(A))])
    R = np.array([[Fp(x) for x in B[k]] for k in range(len(B))])
    O_ = np.array([[Fp(x) for x in C[k]] for k in range(len(C))])
    return L, R, O_, r
    # L = np.transpose(L)


print(testing())

print(L.dot(r) * R.dot(r) - O_.dot(r))

print(r)
print(A)
print(B)
print(C)

"""
def use_Pinocchio_paper_example(b0, b1, b2, b3):
    a = np.array([Fp(x) for x in [1, b0, b1, b2, b3]])
    L0 = np.array([Fp(x) for x in [2, 0, 0, 0, 0]])
    L1 = np.array([Fp(x) for x in [0, 1, 0, 0, 0]])
    L2 = np.array([Fp(x) for x in [0, 0, 1, 0, 0]])
    L3 = np.array([Fp(x) for x in [0, 0, 0, 1, 0]])
    L4 = np.array([Fp(x) for x in [0, 0, 0, 0, 1]])
    L = np.array([L0, L1, L2, L3, L4])

    R0 = np.array([Fp(x) for x in [1, 0, 0, 0, 0]])
    R1 = np.array([Fp(x) for x in [0, 1, 0, 0, 0]])
    R2 = np.array([Fp(x) for x in [0, 0, 1, 0, 0]])
    R3 = np.array([Fp(x) for x in [0, 0, 0, 1, 0]])
    R4 = np.array([Fp(x) for x in [0, 0, 0, 0, 1]])
    R = np.array([R0, R1, R2, R3, R4])

    O0 = np.array([Fp(x) for x in [0, 1, 2, 4, 8]])
    O1 = np.array([Fp(x) for x in [0, 1, 0, 0, 0]])
    O2 = np.array([Fp(x) for x in [0, 0, 1, 0, 0]])
    O3 = np.array([Fp(x) for x in [0, 0, 0, 1, 0]])
    O4 = np.array([Fp(x) for x in [0, 0, 0, 0, 1]])
    O_ = np.array([O0, O1, O2, O3, O4])

    print("La * Ra: = ", L.dot(a) * R.dot(a))
    print("Oa = ", O_.dot(a))
    return L, R, O_, a
"""
