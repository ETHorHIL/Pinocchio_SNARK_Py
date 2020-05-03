import numpy as np

# a = np.array([Fp(x) for x in [1, wire_a, wire_b, wire_NAND]])
U1 = np.array([x for x in [-1, 2, 0, 0]])
U2 = np.array([x for x in [-1, 0, 2, 0]])
U3 = np.array([x for x in [-1, 0, 0, 2]])
U4 = np.array([x for x in [-5, 2, 2, 4]])
L = np.array([U1, U2, U3, U4])
R = L
O_ = L * R
print("L*R: ")
print(O_)
print(O_.dot([1, 0, 1, 1]))
