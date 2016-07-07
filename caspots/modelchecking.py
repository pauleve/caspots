#!/usr/bin/env python

import subprocess

from caspots.dataset import *
from caspots.model import *

U_GENERAL = "general"
U_ASYNC = "asynchronous"

MODES = [U_GENERAL, U_ASYNC]

def make_smv(dataset, bnmodel, destfile, update=U_GENERAL):

    constants = dataset.stimulus.difference(bnmodel.nodes)

    clampable = bnmodel.nodes.intersection(dataset.inhibitors.union(dataset.stimulus))

    vs = set([f.var for f in bnmodel.formula.values()])

    constants.update(bnmodel.nodes.difference(vs))

    smv = open(destfile, "w")
    smv.write("MODULE main\n")
    smv.write("\nVAR\n")
    smv.write("\tstart: boolean;\n")
    for n in constants:
        smv.write("\tn_%s: boolean;\n" % n)
    for n in vs:
        smv.write("\tn_%s: boolean;\n" % n)
        smv.write("\tu_%s: boolean;\n" % n)
        if n in clampable:
            smv.write("\tC_%s: {0,1,-1};\n" % n)

    smv.write("\nASSIGN\n")
    smv.write("next(start) := FALSE;\n")
    for n in constants:
        smv.write("next(n_%s) := n_%s;\n" % (n,n))
    for n in vs:
        smv.write("next(n_%s) := case " % n)
        if n not in dataset.readout:
            smv.write("start: {TRUE, FALSE}; ")
        smv.write("u_%s: F_%s; TRUE: n_%s; esac;\n" % (n, n, n))
        if n in clampable:
            smv.write("next(C_%s) := C_%s;\n" % (n, n))
        #smv.write("next(u_%s) := {TRUE, FALSE};\n" % n)
    smv.write("\nDEFINE\n")
    for f in bnmodel.formula.values():
        n = f.var
        if n in clampable:
            smv.write("F_%s := case C_%s=0: %s; " % (n, n, f.nusmv_expr()))
            smv.write("C_%s=1: TRUE; C_%s=-1: FALSE; esac;\n" % (n, n))
        else:
            smv.write("%s\n" % str(f))

    for exp in dataset.experiments.values():
        setup = []
        for (n, c) in exp.mutations.items():
            if n not in constants and n not in clampable:
                continue
            setup.append("%sn_%s" % ("!" if c < 0 else "", n))
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
                if n not in bnmodel.nodes and n not in dataset.stimulus:
                    continue
                state.append("%sn_%s" % ("!" if not v else "", n))
            smv.write("E%d_T%d := %s;\n" % (exp.id, t, " & ".join(state)))

    fpconds = ["n_%s = F_%s" % (n, n) for n in vs]
    smv.write("FIXEDPOINTS := %s;\n" % " & ".join(fpconds))

    smv.write("\nTRANS\n")
    smv.write("  next(start) != start")
    for n in vs:
        smv.write("\n| next(n_%s) != n_%s" % (n,n))
        smv.write("\n| next(u_%s) != u_%s" % (n,n))
    smv.write("\n| FIXEDPOINTS")
    smv.write(";\n")

    if update == U_ASYNC:
        for n in vs:
            cond = " & ".join(["!u_%s" % m for m in vs if m != n])
            smv.write("TRANS u_%s -> %s;\n" % (n, cond))

    smv.write("\nINIT\n")
    smv.write("(start")
    for n in vs:
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

def verify(dataset, bnmodel, destfile, *args, **kwargs):
    smvfile = make_smv(dataset, bnmodel, destfile, *args, **kwargs)
    output = subprocess.check_output(["NuSMV", "-dcx", smvfile])
    ret = output.strip().split()[-1].decode()
    return ret == "true"

if __name__ == "__main__":
    from argparse import ArgumentParser
    import os
    import tempfile
    parser = ArgumentParser()
    parser.add_argument("dataset")
    parser.add_argument("bnmodels")
    parser.add_argument("--asynchronous", action="store_true", default=False)
    args = parser.parse_args()

    update = U_ASYNC if args.asynchronous else U_GENERAL

    name = os.path.basename(args.dataset).replace(".lp", "")
    dataset = Dataset(name)
    dataset.feed_from_asp(args.dataset)
    for bnmodel in iter_bnmodels(open(args.bnmodels)):
        fd, smvfile = tempfile.mkstemp(".smv")
        os.close(fd)
        #smvfile = "/tmp/bn.smv"
        ret = verify(dataset, bnmodel, smvfile, update=update)
        print("%s/%06d %s" % (name, bnmodel.id, ret))
        os.unlink(smvfile)

