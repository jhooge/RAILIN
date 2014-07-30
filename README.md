RAILIN
======

Resonance Assignment by Integer Linear Programming

## Usage

```sh
usage: railin.py [-h] [--reffile REFFILE] [-a ASSIGNMENTS] [-t TOLERANCE]
               [-s {Joint,JointNew,CplexILP,ILP,Single}] [-v]
               resfile seqfile

positional arguments:
  resfile               Path to list of pseudo-residues in ASCII format.
  seqfile               Path to sequence file in FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  --reffile REFFILE     Path to BMRB reference shifts in ASCII format.
  -a ASSIGNMENTS, --assignments ASSIGNMENTS
                        Number of assignments to be computed (default: 100).
  -t TOLERANCE, --tolerance TOLERANCE
                        Linking Tolerance in ppm (default: 0.6).
  -s {Joint,JointNew,CplexILP,ILP,Single}, --strategy {Joint,JointNew,CplexILP,ILP,Single}
                        Assignment strategy (default: CplexILP).
  -v, --verbose         Turn on verbosity option.
```

## Introduction
RAILIN is a python toolkit for an automated amino acid side chain assignment from nuclear magnetic resonance frequencies. This problem can be formulated as linear assignment problem and solved by constrained optimization techniques. RAILIN uses a MAP (maximum aposteriori) classification approach to estimate the probability, that the resonance frequencies of each spin-system correspond to one of the 21 amino acids. Prior knowledge is provided by the protein primary structure. This amino acid sequence is used to formulate a constrained linear assignment problem, which in turn is solved using Integer Linear Programming. An integer programming problem is a mathematical optimization or feasibility program in which some or all of the variables are restricted to be integers. In many settings the term refers to integer linear programming (ILP), in which the objective function and the constraints (other than the integer constraints) are linear.  

## Installation
### Standard Installation
#### Dependencies
RAILIN depends on the following packages:

- `numpy` -- required ([numpy.scipy.org](http://www.numpy.org/))
- `scipy` -- required ([scipy.org](http://www.scipy.org/))
- `matplotlib` -- required ([matplotlib.org](http://matplotlib.org/))
- `networkx` -- required ([networkx.github.io/](https://networkx.github.io/))
- `IBM ILOG CPLEX Optimizer` -- required ([cplex-optimizer](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/))

To check whether any of the above packages are already installed on your system, try:

```python
$ python
>>> import numpy
>>> import scipy
>>> import matplotlib
>>> import networkx
>>> import cplex
```
ImportError means the package is not installed.


### Install dependencies on unix-like systems

1. Use your package manager to install the required python packages; for example on Debian/Ubuntu:

```sh
$ sudo apt-get install python-numpy
$ sudo apt-get install python-scipy
$ sudo apt-get install python-matplotlib
$ sudo apt-get install python-networkx
```

2. Download and install the IBM ILOG CPLEX Optimizer:

For testing purposes a trial version can be found [here](http://www14.software.ibm.com/webapp/download/search.jsp?pn=IBM+ILOG+CPLEX). For full functionality you need to register at ibm.com and get a licenced copy.
