import sys
from solvers import generate
import argparse
import json

parser = argparse.ArgumentParser()

parser.add_argument("outfile", type=str, help="Fichier d'output")
parser.add_argument("dimension", type=int, help="Dimension du problème à génerer (= nombre de variables)")
parser.add_argument("-i","--inequality", action="store_true",help="Utilise la généralisation Gx <= h")
parser.add_argument("-p","--precision", action="store_true",help="Utilise b sans incertitude")
parser.add_argument("-s","--sample",type=int,help="Ajoute les échantillons bs à la configuration")
args = parser.parse_args()

# les arguments par defaut de (?)
borne = 'bornes'
precision = False
ech = 0

if args.inequality:
    borne = 'inq'

if args.precision:
    precision = True

if args.sample != None:
    ech = args.sample

if args.sample != None and args.inequality:
    A, b, Cb, C, d, G, h, x, bs = generate(args.dimension, borne,precision,ech)
    data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)),
          "C" : list(map(list, C)), "d": list(d), "G" : list(G), "G": list(h), "x" : list(x),"bs" : list(map(list,bs))}
elif args.sample != None:
    A, b, Cb, C, d, l, u, x, bs = generate(args.dimension, borne,precision,ech)
    data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)),
          "C" : list(map(list, C)), "d": list(d), "l" : list(l), "u": list(u), "x" : list(x),"bs" : list(map(list,bs))}
elif args.inequality:
    A, b, Cb, C, d, G, h, x = generate(args.dimension, borne,precision,ech)
    data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)),
          "C" : list(map(list, C)), "d": list(d), "G" : list(G), "G": list(h), "x" : list(x)}
else:
    A, b, Cb, C, d, l, u, x = generate(args.dimension, borne,precision,ech)
    data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)),
          "C" : list(map(list, C)), "d": list(d), "l" : list(l), "u": list(u), "x" : list(x)}

with open(args.outfile, "w") as fd:
    fd.write(json.dumps(data))

print("Problème généré avec succès.")

