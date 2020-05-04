from finitefield.finitefield import FiniteField
from finitefield.polynomial import polynomialsOver
from ssbls120 import Fp, Poly, Group
import numpy as np
import random

G = Group.G
GT = Group.GT

ROOTS = [Fp(i) for i in range(128)]


def random_fp():
    return Fp(random.randint(0, Fp.p-1))


def vanishing_poly(S):
    """
    args:
        S (m vector)
    returns:
        p(X) = (X -S1)*(X-S2)* ... * (X-Sm)

    """
    p = Poly([Fp(1)])
    for s in S:
        p *= Poly([-s, Fp(1)])
    return p


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


def use_babySNARK_paper_example(wire_a, wire_b, wire_NAND):
    a = np.array([Fp(x) for x in [1, wire_a, wire_b, wire_NAND]])
    U1 = np.array([Fp(x) for x in [-1, 2, 0, 0]])
    U2 = np.array([Fp(x) for x in [-1, 0, 2, 0]])
    U3 = np.array([Fp(x) for x in [-1, 0, 0, 2]])
    U4 = np.array([Fp(x) for x in [-5, 2, 2, 4]])
    L = np.array([U1, U2, U3, U4])
    R = L
    O1 = np.array([Fp(x) for x in [1, 0, 0, 0]])
    O_ = np.array([O1, O1, O1, O1])

    print("Using Paper Example...")
    print("a: ", a)
    print("L: ")
    print(L)
    print("R: ")
    print(R)
    print("O: ")
    print(O_)
    (m, n) = L.shape
    print("Matrix shape: ", m, " X ", n)
    print("This assignment is satisfying?",
          (L.dot(a) * R.dot(a) == O_.dot(a)).all())
    print("La: ")
    print(L.dot(a))
    print("Ra: ")
    print(R.dot(a))
    print("La * Ra: ")
    print(L.dot(a) * R.dot(a))
    print("Oa: ")
    print(O_.dot(a))
    return L, R, O_, a


def evaluate_in_exponent(powers_of_tau, poly):
    # powers_of_tau:
    #     [G*0, G*tau, ..., G*(Tau**m)]
    # poly:
    #    degree m-bound polynomial in coefficient form
    print("P.degree: ", poly.degree())
    print("Powers of tau: ", len(powers_of_tau))
    assert poly.degree() + 1 < len(powers_of_tau)
    return sum([powers_of_tau[i] * poly.coefficients[i] for i in
                range(poly.degree()+1)], G*0)


# setup
def non_zk_SNARK_setup(L, R, _O, n_stmt):
    """
    L, R, O: the matrix representing the problem equations
    n_stmt: number of entries
    L,R,O v_one: validators restriction on variables
    """
    (m, n) = L.shape

    assert L.shape == R.shape == _O.shape
    assert n_stmt < n

    # Generate roots for each gate
    # Make sure there are roots for each row (equation)
    # Doesnt matter what values the roots have
    # Roots are public
    global ROOTS
    if len(ROOTS) < m:
        ROOTS = tuple(range(m))

    # Generate polynomials for columns of L, R, O
    # intrerpolate for points (x,y) where x's are the roots and y's are the
    # values of the k-th row
    # This is public
    Ls = [Poly.interpolate(ROOTS[:m], L[:, k]) for k in range(n)]
    Rs = [Poly.interpolate(ROOTS[:m], R[:, k]) for k in range(n)]
    Os = [Poly.interpolate(ROOTS[:m], _O[:, k]) for k in range(n)]

    # Sample Random Trapdoors
    # These are only known to the trusted party that generates the setup
    global s, rho_l, rho_r, rho_o, alpha_l, alpha_r, alpha_o, beta, gamma
    s = random_fp()
    # rho_l and rho_r rho_o for randomization of generators
    rho_l = random_fp()
    rho_r = random_fp()
    # shift parameters for polynomials alpha and beta
    alpha_l = random_fp()
    alpha_r = random_fp()
    alpha_o = random_fp()
    beta = random_fp()
    gamma = random_fp()

    # Set rho_o and the operand generators
    rho_o = rho_l * rho_r
    g_l = G * rho_l
    g_r = G * rho_r
    g_o = G * rho_o

    # Set the proving keys
    powers_of_s = [G * (s ** i) for i in range(m+1)]

    # I am unsure of I should add the l_0 terms or if it is covered by the
    # n_stmt
    # for the example we have n = 3 variables and 4 polynomials which already
    # includes the 1 polynomial
    # l_0 = Poly.interpolate(Fp(0), L_v_one)
    # r_0 = Poly.interpolate(Fp(0), R_v_one)
    # o_0 = Poly.interpolate(Fp(0), O_v_one)
    # g_l_to_l_i = [g_l * l_o] + [g_l * L_i(s) for L_i in Ls]

    gl_to_li = [g_l * L_i(s) for L_i in Ls]
    gr_to_ri = [g_r * R_i(s) for R_i in Rs]
    go_to_oi = [g_o * O_i(s) for O_i in Os]

    gl_to_li_shift_a = [g_l * (L_i(s) * alpha_l) for L_i in Ls[n_stmt:]]
    gr_to_ri_shift_a = [g_r * (R_i(s) * alpha_r) for R_i in Rs[n_stmt:]]
    go_to_oi_shift_a = [g_o * (O_i(s) * alpha_o) for O_i in Os[n_stmt:]]

    gl_to_li_shift_b = [g_l * (L_i(s) * beta) for L_i in Ls[n_stmt:]]
    gr_to_ri_shift_b = [g_r * (R_i(s) * beta) for R_i in Rs[n_stmt:]]
    go_to_oi_shift_b = [g_o * (O_i(s) * beta) for O_i in Os[n_stmt:]]

    # Leaving out the ZK part which would require the encrypted target and shft

    proving_key = [powers_of_s,
                   [gl_to_li, gr_to_ri, go_to_oi],
                   [gl_to_li_shift_a, gr_to_ri_shift_a, go_to_oi_shift_a],
                   [gl_to_li_shift_b, gr_to_ri_shift_b, go_to_oi_shift_b]]

    # verification key
    t = vanishing_poly(ROOTS[:m])
    # this could be a problem: in the other version it was G not g_o below
    go_to_t = g_o * t(s)
    statemet_polys = [gl_to_li[:n_stmt], gr_to_ri[:n_stmt], go_to_oi[:n_stmt]]
    verifier_key = [go_to_t, statemet_polys,
                    G * alpha_l,
                    G * alpha_r,
                    G * alpha_o,
                    G * gamma,
                    G * (gamma * beta)]

    return proving_key, verifier_key, [Ls, Rs, Os]


# Prover
def non_zk_SNARK_prover(L, R, O_, LROpoly, n_stmt, proving_key, a):
    """
    U: the matrix m*n representing the problem equations
    n_stmt: the first l entries of the solution vectore representing the stmt
    CRS: the common reference string, babysnark_setup()[0]
    a: the vector [solution + witness]
    """
    (m, n) = L.shape
    assert L.shape == R.shape == O_.shape
    assert n == len(a)
    print("n_stmt", n_stmt)
    # assert len(CRS) == (m+1) + 2 + (n - n_stmt)
    print("len(proving_key): ", len(proving_key))
    assert len(proving_key) == 4
    assert len(proving_key[0]) == m + 1
    assert len(proving_key[1]) == 3
    assert len(proving_key[1][0]) == len(proving_key[1][1]) \
        == len(proving_key[1][2]) == n
    assert len(proving_key[2]) == 3
    assert len(proving_key[2][0]) == len(proving_key[2][1]) \
        == len(proving_key[2][2]) == n - n_stmt
    assert len(proving_key[3]) == 3
    assert len(proving_key[3][0]) == len(proving_key[3][1]) \
        == len(proving_key[3][2]) == n - n_stmt

    assert len(ROOTS) >= m

    # parse the proving key
    powers_of_s = proving_key[0]

    gl_to_li = proving_key[1][0]
    gr_to_ri = proving_key[1][1]
    go_to_oi = proving_key[1][2]

    gl_to_li_shift_a = proving_key[2][0]
    gr_to_ri_shift_a = proving_key[2][1]
    go_to_oi_shift_a = proving_key[2][2]

    gl_to_li_shift_b = proving_key[3][0]
    gr_to_ri_shift_b = proving_key[3][1]
    go_to_oi_shift_b = proving_key[3][2]

    Ls = LROpoly[0]
    Rs = LROpoly[1]
    Os = LROpoly[2]

    # Target is the vanishing polynomial
    t = vanishing_poly(ROOTS[:m])

    # this is missing l_0, r_0, o_0
    L_big = Poly([])
    for k in range(n):
        L_big += Ls[k] * a[k]

    R_big = Poly([])
    for k in range(n):
        R_big += Rs[k] * a[k]

    O_big = Poly([])
    for k in range(n):
        O_big += Os[k] * a[k]

    # Finally p
    p = L_big * R_big - O_big

    # compute the H term, i.e. cofactor H so that P = T * H
    h = p/t
    assert p == h*t

    # assign provers variable to encrypted polynomials
    # TODo the one at the end might have to be changed s l_0 for all
    print("checking for typeerror: ")
    for k in range(n_stmt, n):
        print("k= ", k)
        print("gl[k]: ", gl_to_li[k])
        print("a[k]:", a[k])
        print("gl * a:", gl_to_li[k] * a[k])

    gLbig_at_s = sum([gl_to_li[k] * a[k] for k in range(n_stmt, n)], G*0)
    gRbig_at_s = sum([gr_to_ri[k] * a[k] for k in range(n_stmt, n)], G*0)
    gObig_at_s = sum([go_to_oi[k] * a[k] for k in range(n_stmt, n)], G*0)

    print("inspection:")
    print("gl_to_li_shift_a lenght: ", len(gl_to_li_shift_a))
    print("a lenght: ", len(a))
    gLbig_at_s_shift = sum([gl_to_li_shift_a[k-n_stmt] * a[k]
                           for k in range(n_stmt, n)], G*0)
    gRbig_at_s_shift = sum([gr_to_ri_shift_a[k-n_stmt] * a[k]
                           for k in range(n_stmt, n)], G*0)
    gObig_at_s_shift = sum([go_to_oi_shift_a[k-n_stmt] * a[k]
                           for k in range(n_stmt, n)], G*0)

    # assign the variable values consistency polynomials
    # from here I have no idea if G * 0 is to be added
    g_to_Z = sum([(gl_to_li_shift_b[k-n_stmt] + gr_to_ri_shift_b[k-n_stmt]
                 + go_to_oi_shift_b[k-n_stmt])
                 * a[k] for k in range(n_stmt, n)], G * 0)

    g_h_at_s = evaluate_in_exponent(powers_of_s, h)

    pi = [gLbig_at_s, gRbig_at_s, gObig_at_s,
          g_h_at_s,
          gLbig_at_s_shift, gRbig_at_s_shift, gObig_at_s_shift,
          g_to_Z]

    return pi


def non_zk_SNARK_verifier(L, R, O_, m, n, verifier_key, a_stmt, pi):
    """
    U: the matrix m*n representing the problem equations
    CRS: the common reference string, babysnark_setup()[0]
    Precomp: precomputation provided by babysnark_setup()[1]
    a_stmt: the first part of the solution vecor, part of the statement
    pi: proof, output of  prover, H, Bw, Vw
    """
    # (m, n) = L[0].shape
    # assert L[0].shape == R[1].shape == O_[2].shape
    # print("     (m,n): ", m, ", ", n)
    # (H, Bw, Vw) = pi

    # parse proof
    gLp = pi[0]
    gRp = pi[1]
    gOp = pi[2]
    gh = pi[3]
    gLshift = pi[4]
    gRshift = pi[5]
    gOshift = pi[6]
    gZ = pi[7]

    # parse verifier_key
    go_to_t = verifier_key[0]
    statemet_polys = verifier_key[1]
    g_alpha_l = verifier_key[2]
    g_alpha_r = verifier_key[3]
    g_alpha_o = verifier_key[4]
    g_gamma = verifier_key[5]
    g_gammabeta = verifier_key[6]

    assert len(ROOTS) >= m
    n_stmt = len(a_stmt)

    gLv = sum([statemet_polys[0][k] * a_stmt[k] for k in range(n_stmt)], G * 0)
    gRv = sum([statemet_polys[1][k] * a_stmt[k] for k in range(n_stmt)], G * 0)
    gOv = sum([statemet_polys[2][k] * a_stmt[k] for k in range(n_stmt)], G * 0)

    # Check 1: variable polynimals consistency
    pair_gLp_g_alpha_l = gLp.pair(g_alpha_l)
    pair_gLshift_g = gLshift.pair(G)

    pair_gRp_g_alpha_r = gRp.pair(g_alpha_r)
    pair_gRshift_g = gRshift.pair(G)

    pair_gOp_g_alpha_o = gOp.pair(g_alpha_o)
    pair_gOshift_g = gOshift.pair(G)

    print("Pairings: gLp_g_alpha_l == gLshift_g: ",
          pair_gLp_g_alpha_l == pair_gLshift_g)
    print("Pairings: pair_gRp_g_alpha_r == pair_gRshift_g: ",
          pair_gRp_g_alpha_r == pair_gRshift_g)
    print("Pairings: pair_gOp_g_alpha_o == pair_gOshift_g: ",
          pair_gOp_g_alpha_o == pair_gOshift_g)

    # Check 2: variable values consistency
    g_P = gLp + gRp + gOp
    pair_g_P_and_g_gammabeta = g_P.pair(g_gammabeta)
    pair_gZ_and_g_gamma = gZ.pair(g_gamma)
    print("Pairins: pair_g_P_and_g_gammabeta == gZ_and_g_gamma: ",
          pair_g_P_and_g_gammabeta == pair_gZ_and_g_gamma)

    # Check 3:  valid operations
    gL = gLp + gLv
    gR = gRp + gRv
    gO = gOp + gOv

    pair_gL_gR = gL.pair(gR)
    pair_go_to_t_gh = go_to_t.pair(gh)
    pair_gO_g = gO.pair(G)
    print("Pairing: pair_gL_gR == pair_go_to_t_gh * pair_gO_g: ",
          pair_gL_gR == pair_go_to_t_gh * pair_gO_g)

    return True


def testingProof():

    n_stmt = 1
    # L, R, O_, a = use_paper_example(1, 1, 0)
    L, R, O_, a = use_Pinocchio_paper_example(0, 1, 0, 0)
    (m, n) = L.shape
    a_stmt = a[:n_stmt]
    print("a_stmt: ", a_stmt)
    print("MxN: ", m, "*",  n, "=", m * n)

    # setup
    print("Computing setup...")
    proving_key, verifier_key, LROpoly = non_zk_SNARK_setup(L, R, O_, n_stmt)
    print("proving_key lenght: ", len(proving_key))
    print("proving_key: ", proving_key)
    print("verifier_key lenght: ", len(verifier_key))
    print("verifier_key: ", verifier_key)
    print("[Ls, Rs, Os]: ", len(LROpoly))
    print("[Ls, Rs, Os]: ", LROpoly)

    # prover
    print("Proving...")
    pi = non_zk_SNARK_prover(L, R, O_, LROpoly, n_stmt, proving_key, a)
    print("pi lenght: ", len(pi))

    print("Verifying...")
    success = non_zk_SNARK_verifier(L, R, O_, m, n,
                                    verifier_key, a_stmt[:n_stmt], pi)
    print("VERIFICATION SUCCESSFUL?", success)


testingProof()
