#!/usr/bin/env/python3

import numpy as np
import scipy as sp
import cvxpy as cp
import sympy as sym
import copy
import tqdm

np.set_printoptions(precision=3, linewidth=180)

def generate(m, borne='bornes'):
    """
    On génère les matrice correspondant à un pb à m variables
    systeme = 0 --> n = m -k , systeme = 1 n = m//2 -1
    """

    n = m - 2
    A = np.random.randn(n, m)
    C = np.random.randn(n, m)
    x = np.random.randn(m)
    b = A @ x
    d = C @ x
    vb = np.array([0.1, 0.05, 0.06])**2 # vecteur des variances de chaque terme du vecteur b
    # vb = (np.random.sample(n)/10)**2

    Cb = np.diag(vb) # matrice de covariance de b
    b_soumis = b + np.sqrt(vb) * np.random.randn(n) # b donné, une réalisation de b
    print(x)

    if borne == 'bornes':
        l = x - np.ones(m) * 50
        u = x + np.ones(m) * 50
        return A, b, C, d, l, u
    else:
        G = np.vstack((np.identity(m), -np.identity(m)))
        h = G @ x + np.random.randn(2 * m)
    return A, b_soumis, Cb, C, d, G, h

def test_sous_determine(m, A, C):
    rank = np.linalg.matrix_rank(A) + np.linalg.matrix_rank(C)
    return rank < m

def systeme_sous_determine(m,A,b,C,d,l,u):

    xl = sym.symbols(list('x{} '.format(i) for i in range(m)))
    name = {}
    for i in range(m):
        name[xl[i]] = i

    AC = np.vstack((A, C))
    bd = np.hstack((b, d))

    x_sol = sym.solve(AC @ xl - bd, xl)
    for k,v in x_sol.items():
        xl[name[k]] = v

    for i in range(m):
        if isinstance(xl[i], sym.numbers.Float):
            l[i] = xl[i]
            u[i] = xl[i]

    for i in range(m):
        if isinstance(xl[i], sym.add.Add):
            a = xl[i].as_coefficients_dict()
            al = (list(a.keys()))
            val_l = []
            val_u = []
        for k in al:
            if k != 1:
                val_l.append(l[name[k]])
                val_u.append(u[name[k]])
            else:
                val_l.append(1)
                val_u.append(1)

        b = np.vstack((al, val_l)).T
        c = np.vstack((al, val_u)).T

        new_l = xl[i].subs(b)
        new_u = xl[i].subs(c)
        if new_l > l[i]:
            l[i] = new_l
        if new_u < u[i]:
            u[i] = new_u

    print("Voici l'intervalle de solutions estimé :")
    for i in range(m):
        print('{} <= {} <= {}'.format(l[i], xl[i], u[i]))

    print("Voici une solution approchée parmi les possibles")
    xc = cp.Variable(m)

    # par méthode de cvxpy
    constraints = [C @ xc == d, l <= xc, xc <= u]
    objective = cp.Minimize(cp.sum_squares(A @ xc - b))
    prob = cp.Problem(objective, constraints)

    # resolution
    prob.solve()
    print("x_cvxpy = {}".format(xc.value))

def solver_lagrange_simple_b(m, A, b, C, d, l, u):

    xl = sym.symbols(list('x{} '.format(i) for i in range(m)))

    objective = np.transpose(A @ xl - b) @ (A @ xl - b)
    constraints = C @ xl - d
    lam = sym.symbols(list('lambda{}'.format(i) for i in range(len(constraints))))
    constraints = lam @ constraints

    lagrange = sym.expand(objective - constraints)
    var = xl + lam

    gradient = []

    # On construit le vecteur gradient comme dérivé selon chaque variable
    # Puis les équations associées
    for i in range(len(var)):
        gradient.append(sym.diff(lagrange, var[i]))

    # On cherche la solution de grandient = 0, cela donne un dictionnaire avec les valeurs pour chaque xl
    x_lagrange = list(sym.solve(gradient, var).values())
    l1 = x_lagrange[:m] >= l
    u1 = x_lagrange[:m] <= u
    if False in l1 or False in u1:
        print("Les bornes ne sont pas entièrement respectées")
    print(x_lagrange[:m])

def solver_lagrange_simple_i(m, A, b, C, d, G, h):
    xl = sym.symbols(list('x{} '.format(i) for i in range(m)))

    objective = np.transpose(A @ xl - b) @ (A @ xl - b)
    constraints = C @ xl - d
    lam = sym.symbols(list('lambda{}'.format(i) for i in range(len(constraints))))
    constraints = lam @ constraints

    lagrange = sym.expand(objective - constraints)
    var = xl + lam

    gradient = []

    # On construit le vecteur gradient comme dérivé selon chaque variable
    # Puis les équations associées
    for i in range(len(var)):
        gradient.append(sym.diff(lagrange, var[i]))

    # On cherche la solution de grandient = 0, cela donne un dictionnaire avec les valeurs pour chaque xl
    x_lagrange = list(sym.solve(gradient, var).values())
    test = G @ x_lagrange[:m] <= h

    if False in test:
        print("Les bornes ne sont pas entièrement respectées")
    print(x_lagrange[:m])

def LuToGh(l, u):
    G = np.vstack((np.identity(m), -np.identity(m)))
    h = np.vstack((u, -l))
    return G,h

def solver_cvxpy(m, A, b, C, d, G, h):
    xc = cp.Variable(m)
    if G.shape != (m, 2 * m):
        G,h = LuToGh(G, h)

    constraints = [C @ xc == d, G @ xc <= h]
    objective = cp.Minimize(cp.sum_squares(A @ xc - b))
    prob = cp.Problem(objective, constraints)
    prob.solve()
    print("x_cvxpy = {}".format(xc.value))



m = 5
A, b, C, d, l, u = generate(m)
solver_lagrange_simple_b(m, A, b, C, d, l, u)
