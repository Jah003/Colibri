iport argparse
from solvers import *

parser = argparse.ArgumentParser()
parser.add_argument('configtype', help="Type de configuration (Chargement depuis un fichier ou génération aléatoire", choices = ["load", "generate"])
parser.add_argument('-d', '--dimension', type=int, help="Dimension du problème à générer")
parser.add_argument('-c', '--configfile', type=str, help="Fichier de configuration à charger")

parser.add_argument('solver', type=str, help="Solver a utiliser", choices=["lagrange_simple_i", "lagrange_simple_b","partition"])
parser.add_argument('-n', '--npart', type=int, help="Nombre de partitions", default=2, required=False)

args = parser.parse_args()

# Partie récupération de la configuration voulu (génération aléatoire ou fichier de config)
if args.configtype == "generate":
    if not args.dimension:
        print("Une dimension doit être donner pour générer un problème aléatoire (-d <dim>)")
        sys.exit(1)

    m = args.dimension
    A, b, C, d, l, u = generate(m)

elif args.configtype == "load":
    if not args.configfile:
        print("Une fichier de configuration doit être donné pour être chargé (-c <fichier de conf>)")
        sys.exit(1)

    data = load_config(args.configfile)

    if not data:
        print(f"La configuration du fichier {args.configfile} à été chargée avec succès.")
    else:
        print("Une erreur est survenue en chargeant la configuration.")
        sys.exit(1)

# Parti appel du solver
required_keys = ["A","b","C","d","l","u"]
err = False
for key in required_keys:
    if not key in data.keys():
        print("Il manque {} dans le fichier de config ou à la génération.".format(key))
        err = True

if not err:
    A = data["A"]
    b = data["b"]
    C = data["C"]
    d = data["d"]
    l = data["l"]
    u  = data["u"]

    if args.solver == "lagrange_simple_b":
            x = solver_lagrange_simple_b(m, A, b, C, d, l, u)

    elif args.solver == "lagrange_simple_i":
        G,h = LuToGh(l,u)
        x = solver_lagrange_simple_i(m, A, b, C, d, G, h)

    elif args.solver == "partition":
        x = solver_partitions(m, A, b, C, d, l, u, args.npart)

