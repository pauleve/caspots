#!/usr/bin/env python

import subprocess

U_GENERAL = "general"
U_ASYNC = "asynchronous"

MODES = [U_GENERAL, U_ASYNC]

def make_smv(dataset, network, destfile, update=U_GENERAL):

    # nodes for which a function is defined
    varying_nodes = set([node for node, _ in network.formulas_iter()])

    # nodes with no function (i.e., constant value)
    constants = set(network.variables()).difference(varying_nodes)
    #constants = dataset.stimulus.difference(varying_nodes)

    clampable = varying_nodes.intersection(dataset.inhibitors.union(dataset.stimulus))

    smv = open(destfile, "w")
    smv.write("MODULE main\n")
    smv.write("\nVAR\n")
    smv.write("\tstart: boolean;\n")
    for n in constants:
        smv.write("\tn_%s: boolean;\n" % n)
    for n in varying_nodes:
        smv.write("\tn_%s: boolean;\n" % n)
        smv.write("\tu_%s: boolean;\n" % n)
        if n in clampable:
            smv.write("\tC_%s: {0,1,-1};\n" % n)

    smv.write("\nASSIGN\n")
    smv.write("next(start) := FALSE;\n")
    for n in constants:
        smv.write("next(n_%s) := n_%s;\n" % (n,n))
    for n in varying_nodes:
        smv.write("next(n_%s) := case " % n)
        if n not in dataset.readout:
            smv.write("start: {TRUE, FALSE}; ")
        smv.write("u_%s: F_%s; TRUE: n_%s; esac;\n" % (n, n, n))
        if n in clampable:
            smv.write("next(C_%s) := C_%s;\n" % (n, n))
        #smv.write("next(u_%s) := {TRUE, FALSE};\n" % n)
    smv.write("\nDEFINE\n")

    def nusmv_of_literal((var, sign)):
        return "%sn_%s" % ("!" if sign == -1 else "", var)

    def nusmv_of_clause(clause):
        expr = " & ".join(map(nusmv_of_literal, clause))
        if len(clause) > 1:
            return "(%s)" % expr
        return expr

    def nusmv_of_clauses(clauses):
        if len(clauses) == 0:
            return "FALSE"
        return " | ".join(map(nusmv_of_clause, clauses))

    for n, clauses in network.formulas_iter():
        expr = nusmv_of_clauses(clauses)
        if n in clampable:
            smv.write("F_%s := case C_%s=0: %s; " % (n, n, expr))
            smv.write("C_%s=1: TRUE; C_%s=-1: FALSE; esac;\n" % (n, n))
        else:
            smv.write("F_%s := %s;\n" % (n, expr))

    for exp in dataset.experiments.values():
        setup = []
        # enforce initial state of clamped nodes
        for (n, c) in exp.mutations.items():
            setup.append("%sn_%s" % ("!" if c < 0 else "", n))
        # specify clamping setting
        for n in clampable:
            if n in exp.mutations:
                c = exp.mutations[n]
                setup.append("C_%s=%s" % (n,c))
            else:
                setup.append("C_%s=0" % n)

        smv.write("E%d_SETUP := %s;\n" % (exp.id, " & ".join(setup) or "TRUE"))
        for t, values in exp.obs.items():
            state = []
            for n, v in values.items():
                state.append("%sn_%s" % ("!" if not v else "", n))
            smv.write("E%d_T%d := %s;\n" % (exp.id, t, " & ".join(state)))

    fpconds = ["n_%s = F_%s" % (n, n) for n in varying_nodes]
    smv.write("FIXEDPOINTS := %s;\n" % " & ".join(fpconds))

    smv.write("\nTRANS\n")
    smv.write("  next(start) != start")
    for n in varying_nodes:
        smv.write("\n| next(n_%s) != n_%s" % (n,n))
        smv.write("\n| next(u_%s) != u_%s" % (n,n))
    smv.write("\n| FIXEDPOINTS")
    smv.write(";\n")

    if update == U_ASYNC:
        for n in varying_nodes:
            cond = " & ".join(["!u_%s" % m for m in varying_nodes if m != n])
            smv.write("TRANS u_%s -> %s;\n" % (n, cond))

    smv.write("\nINIT\n")
    smv.write("(start")
    for n in varying_nodes:
        smv.write(" & !u_%s" % n)
    smv.write(");\n")

    def ctl_of_exp(exp):
        ts = list(exp.obs.keys())
        t0 = ts.pop(0)
        ctl = "(E%d_SETUP & E%d_T%d) -> " % (exp.id, exp.id, t0)
        for t in  ts:
            ctl += "EF (E%d_T%d & " % (exp.id, t)
        ctl = ctl[:-2] + ")"*len(ts)
        return "(%s)" % ctl

    smv.write("\nSPEC (\n  ")
    smv.write("\n& ".join([ctl_of_exp(exp) for exp in dataset.experiments.values()]))
    smv.write("\n);\n")
    smv.write("\n")
    smv.close()
    return destfile

def verify(dataset, network, destfile, *args, **kwargs):
    smvfile = make_smv(dataset, network, destfile, *args, **kwargs)
    output = subprocess.check_output(["NuSMV", "-dcx", smvfile])
    ret = output.strip().split()[-1].decode()
    return ret == "true"

