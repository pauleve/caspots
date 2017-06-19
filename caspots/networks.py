
import itertools as it

import gringo

import pandas as pd

from caspo.core import Clause

from .asputils import *

def domain_of_networks(networks, hypergraph, dataset):
    fs = funset(networks)
    domain = ["1{%s}1." % ("; ".join(["model(%d)" % i for i in range(len(networks))]))]

    formulas = set()
    for network in networks:
        formulas = formulas.union(it.imap(lambda (_, f): f, network.formulas_iter()))
    formulas = pd.Series(list(formulas))

    for i, network in enumerate(networks):
        for v, f in network.formulas_iter():
            f= gringo.Fun("formula", [v, formulas[formulas == f].index[0]])
            domain.append("%s :- model(%d)." % (f,i))

    return "%s%s\n" % (fs.to_str(), "\n".join(domain))

def restrict_with_partial_bn(hypergraph, partial_bn_file):
    asp = []
    with open(partial_bn_file) as f:
        for line in f:
            line = line.strip().replace(' ','')
            if not line:
                continue
            a, clauses = line.split("=")
            a = a.strip()
            ai = hypergraph.nodes[hypergraph.nodes == a].index[0]
            clauses = set([Clause.from_str(c) for c in clauses.split("|")])
            for hi in hypergraph.hyper[hypergraph.hyper == ai].index:
                if hypergraph.clauses[hi] in clauses:
                    asp.append("dnf({0},{1}).".format(ai, hi))
                else:
                    asp.append(":- dnf({0},{1}).".format(ai,hi))
    return "\n".join(asp)
