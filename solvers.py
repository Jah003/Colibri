#!/usr/bin/env/python3

import numpy as np
import scipy as sp
import cvxpy as cp
import sympy as sym
import copy
import tqdm
import sys
import json
import argparse
from random import randint

np.set_printoptions(precision=3, linewidth=180)

def load_config(config_file, keys=['A','b','Cb','C','l','u']):
    """
    Parse le fichier de config et retourne les valeurs sérialiser
    Echoue ( et renvoie None) si le fichier n'existe pas, ne contient pas du json valide, ou qu'il manque une ou plusieurs des valeurs.
    """

    try:
        fd = open(config_file, "r")
    except FileNotFoundError:
        print(f"Erreur, le fichier {config_file} n'as pas peu être ouvert en lecture.")
        return None
    try:
        json_data = json.load(fd)
    except json.decoder.JSONDecodeError as e:
        print(f"Erreur, le fichier {config_file} contient du JSON invalide : {e}")
        return None

    values = []
    err = False

    for key in keys:
        if not key in json_data:

            print(f"Erreur, le fichier de config ne contient pas {key}.")

            # On ne return pas ici, cela permet d'affiché la totalité des clés manquantes
            err = True
        else:

            val = np.array(json_data[key])
            values.append(val)

    # On échoue si une ou plusieurs clés sont manquantes
    if err:
        return None

    return values


def generate(m, borne='bornes',precision=False,ech=0):
    """
    m la dimension de x
    nombre de lignes de A et C détermine aleatoirement
    precision : si F retourne un echantillon de b, si T retrouve b_bar
    ech : si > 0 retourne le tableau bs pour iteration
    """

    n = randint(1,m)
    A = np.random.randn(n, m)
    C = np.random.randn((m-n)*randint(1,m), m)
    x = np.random.randn(m)
    b = A @ x
    d = C @ x
    vb = (np.random.sample(n/5)**2 # vecteur des variances de chaque terme du vecteur b

    Cb = np.diag(vb) # matrice de covariance de b
    b_soumis = b + np.sqrt(vb) * np.random.randn(n) # b donné, une réalisation de b
    if precision==True:
        b_soumis = b
    if ech>0:
        bs = b[:, None]+np.sqrt(vb)[:, None]*np.random.randn(n, ech)
    if borne == 'bornes':
        l = x - np.ones(m) * 50
        u = x + np.ones(m) * 50
        if ech>0 :
            return A, b, Cb, C, d, l, u, bs
        else:
            return A, b, Cb, C, d, l, u
    else:
        G = np.vstack((np.identity(m), -np.identity(m)))
        h = G @ x + np.random.randn(2 * m)
        if ech>0 :
            return A, b_soumis, Cb, C, d, G, h, bs
        else:
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
    return x_lagrange[:m]

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
    return x_lagrange[:m]

def LuToGh(l, u):
    m = len(l)
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

def test_diag(a):
  """
  utilise pour construire la matrice var-cov
  T si tous les elements de la diag de a sont entre 0 et 1
  """
  for i in range(a.shape[0]):
    if a[i,i] <= 0 or a[i,i]>=1:
      return False
  return True

def simplifie_sd(m,C,d,l,u,k,x,sd):
  """
  reduit l ecart-type sd pour correspondre au contrainte Cx=d et l<x<u
  x est une solution trouvee
  k est la precision, vaut 1,2, ou 3 (du + precis ou - precis)
  """
  xl = sym.symbols(list('x{} '.format(i) for i in range(m)))
  name = {}
  for i in range(m):
    name[xl[i]] = i
  x_sol = sym.solve(C@xl - d, xl)

  rg = np.linalg.matrix_rank(C)

  sub1 = []
  sub2 = []
  for i in range(rg,m):
    if x[i]-k*sd[i] < l[i]:
      if x[i]+k*sd[i] > u[i]:
        sd[i] = min(np.abs((l[i]-x[i])/k),np.abs((u[i]-x[i])/k))
      else :
        sd[i] = np.abs((l[i]-x[i])/k)
    elif x[i]+k*sd[i] > u[i]:
      sd[i] = np.abs((u[i]-x[i])/k)
    sub1.append([xl[i],x[i]- k*sd[i]])
    sub2.append([xl[i],x[i]+ k*sd[i]])

  for key,v in x_sol.items():
    #print(x[name[key]])
    i = name[key]
    a = v.subs(sub1)
    b = v.subs(sub2)

    if a < l[i]:
      if b > u[i]:
        sd[i] = min(np.abs((l[i]-a)/k),np.abs((u[i]-b)/k))
      else :
        sd[i] = np.abs((l[i]-a)/k)
    elif b > u[i]:
      sd[j] = np.abs((u[i]-b)/k)

  return sd**2


def make_sd(m,A,Cb,C,d,l,u,k,x):
  """
  calcule sd de x pour Ax=b a partir de Cb
  puis retourne le resultat simplifie par methode precedente
  """
  pinvA = np.linalg.pinv(A)
  pinvAt = np.linalg.pinv(A.T)
  prod = pinvA@Cb@pinvAt
  sommeA = np.eye(m) - pinvA@A
  kxa = -1*np.eye(m,m)

  while test_diag(kxa)==False :
    kxa = prod + sommeA@np.random.randn(m,m)
  sq = []
  for i in range(m):
    sq.append(np.sqrt(kxa[i,i]))
  sq = np.array(sq)

  return simplifie_sd(m,C,d,l,u,k,x,sq)

def compare_sd(A,B):
  res = []
  if A.shape==B.shape:
    for i in range(A.shape[0]):
      res.append(np.abs(np.sqrt(A[i,i])-np.sqrt(B[i,i])))
    return res
  else:
    for i in range(A.shape[0]):
      res.append(np.abs(np.sqrt(A[i,i])-np.sqrt(B[i])))
    return res


def creation_partitions(ls, us, npart, dim):
    """
    Cette fonction va nous servir à créer nos différentes partitions composées de m "bebés" intervalles,
    Pour créer toutes les partitions, on crée toutes les combinaisons possibles de "bébés" intervalles
    """

    assert(dim == len(ls) == len(us))

    possible_couples = []

    for m in range(dim):
        possible_couples.append([(li,ui) for li,ui in create(ls[m],us[m],npart)])

    ls_us = itertools.product(*possible_couples) #On produit les différentes combinaisons et on les retourne

    return [list(map(list, zip(*x))) for x in ls_us]

def resolution_base(A,b,C,d,G,h): # Cette fonction resout le probleme sans utiliser la partition

    xes= []

    for b in tqdm.tqdm(bs.T):

        # definition du problème
        xe = cp.Variable(m)

        objective = cp.Minimize(cp.sum_squares((A @ xe) - b))
        constraints = [(C @ xe) == d, G @ xe <= h]
        prob = cp.Problem(objective, constraints)
        # resolution
        try :
            prob.solve()
        except Exception as e :
            print(e)
        if xe.value is None:
            continue
        xes.append(xe.value)

    # print("Dcp ? {}".format(prob.is_dcp()))

    return xes

def solver_partitions(m, A, b, C, d, G, h, n_part, affichage=False, affichage_no_sol=False):
    """
    Fonction principale qui va résoudre le système pour chaque b et chaque partitions
    """

    u = np.array(h[:m])
    lh = np.array(h[-m:])
    hp = []

    if (affichage):

        for li,ui in creation_partitions(lh, u, n_part,m):
            li = np.array(li)
            ui = np.array(ui)
            print("l={},u={}".format(-li,ui))


    min_local = []

    for li,ui in creation_partitions(lh,u,n_part,m):

        li = np.array(li)
        ui = np.array(ui)

        hp = np.hstack([ui,li])
        xe = cp.Variable(m)
        obj = cp.Minimize(cp.sum_squares( (A @ xe) - bi) )
        cstr = [(C @ xe) == d, G @ xe <= hp]
        prob = cp.Problem(obj, cstr)

        try :
            prob.solve()
        except Exception as e :
            print(e)

        # Pas de solution
        if xe.value is None:
            if (affichage_no_sol) :
                print("Pas de solution :( pour {} <= x <= {}".format(-li,ui))
            continue


        # Stock les solutions de chaque partitions
        min_local.append((xe.value, np.linalg.norm( (A @ xe.value) - bi)**2 ))

    # Récupère la solution dont la norme de A@xsol-b^2 est la plus petite
    return min(min_local, key=lambda couple: couple[1]) )

def solver_mean_partitions(m, A, bs, C, d, G, h, n_part):

    """
    Pour chaque b, résout le problème avec le solver par partition
    et renvoie la moyenne des xe obtenues
    """

    xes = []
    for b in bs:
    xe = solve_partition(m, A, b, C, d, G, h, n_part)
        xes.append(xe)

    x_opt1 = xes1.T.mean(axis=1)
    return x_opt1, xes1
