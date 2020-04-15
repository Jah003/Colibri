import sys
from solvers import generate
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("outfile", type=str, help="Fichier d'output")
parser.add_argument("dimension", type=int, help="Dimension du problème (=nombre de variables)")
parser.add_argument("-b", "--bornes", action="store_true")
args = parser.parse_args()

if args.bornes:
    A, b, Cb, C, d, G, h = generate(args.dimension)
else:
    A, b, Cb, C, d, G, h = generate(args.dimension, borne='')

data = {"A" : list(map(list, A)), "b" : list(b), "Cb": list(map(list, Cb)), "C" : list(map(list, C)), "d": list(d), "G" : list(map(list, G)), "h": list(h)}

with open(args.outfile, "w") as fd:
    fd.write(str(data))


