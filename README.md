# Colibri

## Installation

### Dépendances

`pip install scipy numpy sympy cvxpy tqdm`

Nous utilisons les librairies suivantes:

* [CVXPY](https://www.cvxpy.org/)
* [Simpy](https://www.sympy.org/en/index.html)
* [Scipy](https://www.scipy.org/)
* [Numpy](https://numpy.org/)

### Logiciel

Après avoir installé les dépendances, il suffit de `git clone https://github.com/Antoxyde/Colibri`.

## Utilisation

Le programme est composé de 2 script executables, `genconf.py` et `colibri.py`.

### Génération d'une configuration aléatoire avec `genconf.py`.
Il faut lui spécifier une fichier de sortie, et la dimension du problème.

Par exemple , `python genconf.py test.json 5`, génére une configuration de dimension 5 et l'écrit dans le fichier `test.json`.

### Lancement des solvers avec `colibry.py`

TODO
