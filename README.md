# Colibri

## Installation

### Dépendances

`pip install scipy numpy sympy cvxpy`

Nous utilisons les librairies suivantes:

* [CVXPY](https://www.cvxpy.org/)
* [Simpy](https://www.sympy.org/en/index.html)
* [Scipy](https://www.scipy.org/)
* [Numpy](https://numpy.org/)

### Logiciel

Après avoir installer les dépendances, il suffit de `git clone https://github.com/Antoxyde/Colibri`.

## Utilisation

Génération d'une configuration aléatoire avec `genconf.py`.

Il faut lui spécifier une fichier de sortie, et la dimension du problème.

par exemple :

`python genconf.py test.json 5`

## Lancement du solver

