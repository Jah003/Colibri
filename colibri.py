import argparse
from solvers import *

parser = argparse.ArgumentParser()
parser.add_argument('configtype', help="Type de configuration (Chargement depuis un fichier ou génération aléatoire", choices = ["load", "generate"])
parser.add_argument('-d', '--dimension', type=int, help="Dimension du problème à générer")
parser.add_argument('-c', '--configfile', type=str, help="Fichier de configuration à charger")

args = parser.parse_args()

# Partie récupération de la configuration voulu (génération aléatoire ou fichier de config)
if args.configtype == "generate":
    if not args.dimension:
        print("Une dimension doit être donner pour générer un problème aléatoire (-d <dim>)")
        sys.exit(1)

    m = args.dimension
    A, b, C, d, G, h = generate(m)

elif args.configtype == "load":
    if not args.configfile:
        print("Une fichier de configuration doit être donné pour être chargé (-c <fichier de conf>)")
        sys.exit(1)

    vals = load_config(args.configfile)

    if vals:
        A, b, C, d, G, h = vals
        m = len(A)
        print(f"La configuration du fichier {args.configfile} à été chargée avec succès.")
    else:
        print("Une erreur est survenue en chargeant la configuration.")
        sys.exit(1)


# Partie résolution du problème

print("Parti solve")
print(m, A, b, C, d, G, h)

#solver_lagrange_simple_b(m, A, b, C, d, l, u)
