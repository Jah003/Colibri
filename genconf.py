import sys
from solvers import generate
import argparse
import json

parser = argparse.ArgumentParser()

parser.add_argument("outfile", type=str, help="Fichier d'output")
parser.add_argument("dimension", type=int, help="Dimension du problème à génerer (= nombre de variables)")
args = parser.parse_args()

A, b, Cb, C, d, l, u = generate(args.dimension, borne="bornes")
data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)), "C" : list(map(list, C)), "d": list(d), "l" : list(l), "u": list(u)}

with open(args.outfile, "w") as fd:
    fd.write(json.dumps(data))

print("Problème généré avec succès.")

