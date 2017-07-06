
import itertools as it

import gringo

import pandas as pd

from caspo.core import Clause

from .asputils import *

def domain_of_networks(networks):

    hg = networks.hg

    domain = ["1{%s}1." % ("; ".join(["model(%d)" % i for i in range(len(networks))]))]

    for i,network in enumerate(networks):
        m = gringo.Fun("model", [i])
        for v,formula in network.formulas_iter():
            variable = hg.nodes[hg.nodes == v].index[0]
            f = gringo.Fun("formula", [v, variable])
            domain.append("%s :- model(%d)." % (f,i))
            for clause in formula:
                clause_idx = hg.clauses_idx[clause]
                d = gringo.Fun("dnf",[variable, clause_idx])
                domain.append("%s :- %s." % (d,m))
                for variable, sign in clause:
                    c = gringo.Fun("clause", [clause_idx, variable, sign])
                    domain.append("%s :- %s." % (c,m))

    return "%s\n" % "\n".join(domain)

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
