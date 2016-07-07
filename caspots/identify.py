
from __future__ import print_function

import math
import os
from subprocess import *
import sys
import tempfile
import time

try:
    import gringo
except ImportError:
    print("!!! gringo python module is not installed. Some features will fail.",
            file=sys.stderr)

from pyzcasp import asp

from caspots.config import *
from caspots.model import *
from caspots import asputils
from caspots.utils import *

def parse_answer(line):
    return asputils.re_answer.findall(line)

def crunch_data(answer, predicate, factor):
    factor = float(factor)
    data = {
        "obs": {},
        "bin": {}
    }
    keys = set()
    for (p, args) in answer:
        if p in ["obs", predicate]:
            args = asputils.parse_args(args)
            key = tuple(args[:3])
            val = args[3]
            if p == "obs":
                val /= factor
            t = "obs" if p == "obs" else "bin"
            data[t][key] = val
            keys.add(key)
    return (keys, data)

def MSE(cd):
    cum = 0
    keys, data = cd
    for key in keys:
        if key not in data["obs"]:
            continue
        cum += (data["obs"][key] - data["bin"][key])**2
    return math.sqrt(cum)

def count_predicate(answer, predicate):
    return len([p for (p,a) in answer if p == predicate])

class ASPSample:
    def __init__(self, result, opts):
        self.opts = opts
        mode = 0
        for line in result:
            line = line.strip()
            if mode == 0:
                if line.startswith("Answer: "):
                    self.id = int(line[7:])
                    mode = 1
            elif mode == 1:
                self.answer = parse_answer(line)
                mode = 2
            elif mode == 2:
                if line.startswith("Optimization: "):
                    self.optimization = line.split(" ")[1]

    def weight(self):
        return self.optimization

    def asp_exclusion(self):
        predicates = ["formula", "dnf", "clause"]
        if self.opts.enum_traces:
            predicates += ["guessed"]
        clauses = ["%s(%s)" % (p,a) for (p,a) in self.answer if p in predicates]
        if self.opts.family == "all":
            nb_formula = count_predicate(answer, "formula")
            nb_dnf = count_predicate(answer, "dnf")
            nb_clause = count_predicate(answer, "clause")
            clauses += [
                "%d{formula(V,I): node(V,I)}%d" % (nb_formula, nb_formula),
                "%d{dnf(I,J): hyper(I,J,N)}%d" % (nb_dnf, nb_dnf),
                "%d{clause(J,V,B): edge(J,V,B)}%d" % (nb_clause, nb_clause)
            ]
        return ":- %s." % ", ".join(clauses)

    def mse(self):
        cd_measured = crunch_data(self.answer, "measured", self.opts.factor)
        cd_guessed = crunch_data(self.answer, "guessed", self.opts.factor)
        mse0 = MSE(cd_measured)
        mse = MSE(cd_guessed)
        return (mse0, mse)

    def model(self):
        bnmodel = BNModel(self.id)
        bnmodel.feed_from_asp(self.answer)
        return bnmodel

    def trace(self, dataset):
        # rewrite dataset using guessed predicate
        for (p, args) in self.answer:
            if p == "guessed":
                eid, t, node, value = asputils.parse_args(args)
                if node not in dataset.readout:
                    continue
                if dataset.experiments[eid].obs[t][node] != value:
                    #print(((eid,t,node),dataset.experiments[eid].obs[t][node], value), file=sys.stderr)
                    dataset.experiments[eid].obs[t][node] = value
        return dataset

    def answerset(self):
        return asp.AnswerSet(["%s(%s)" % a for a in self.answer])


class ASPSolver:
    def __init__(self, termset, opts):
        self.termset = termset
        self.data = termset.to_str()
        self.opts = opts
        self.debug = opts.debug
        self.domain = [aspf("guessBN.lp")]
        if opts.fully_controllable:
            self.domain.append(aspf("guessBN-controllable.lp"))

    def sample(self, first, *args):
        gringo_cmd = ["gringo"] + self.domain + [
            aspf("supportConsistency.lp"),
            aspf("normalize.lp"),
            aspf("showMeasured.lp")]
        gringo_cmd += ["-"]
        if self.opts.family == "subset":
            gringo_cmd.append(aspf("minimizeSizeOnly.lp"))
        if first:
            gringo_cmd += [aspf("minimizeWeightOnly.lp")]
        gringo_cmd += list(args)
        clasp_cmd = ["clasp", "--conf=trendy", "--stats", "--opt-strat=usc",
                        "--quiet=1"]
        if self.debug:
            dbg("# %s" % (" ".join(gringo_cmd + ["|"] + clasp_cmd)))
        gringo = Popen(gringo_cmd, stdin=PIPE, stdout=PIPE)
        clasp = Popen(clasp_cmd, stdin=gringo.stdout, stdout=PIPE)
        gringo.stdin.write(self.data)
        gringo.stdin.close()
        sample = ASPSample(clasp.stdout, opts=self.opts)
        clasp.stdout.close()
        clasp.wait()
        self.retcode = clasp.returncode
        return sample

    def solution_samples(self):
        i = 1
        if self.debug:
            dbg("# model %d" %i)
        s = self.sample(True)
        yield s

        weight = s.weight()
        fd, excludelp = tempfile.mkstemp(".lp")
        os.close(fd)

        with open(excludelp, "w") as f:
            f.write("%s\n" % s.asp_exclusion())

        args = [aspf("tolerance.lp"), excludelp,
            "-c", "minWeight=%s" % weight,
            "-c", "maxWeight=%s" % weight]

        while True:
            s = self.sample(False, *args)
            if self.retcode in [10, 30]:
                i += 1
                if self.debug:
                    dbg("# model %d" %i)
                yield s
                with open(excludelp, "a") as f:
                    f.write("%s\n" % s.asp_exclusion())
            elif self.retcode == 20:
                print("# Enumeration complete")
                break
            else:
                print("# ERROR")
                break

        os.unlink(excludelp)

    @asp.cleanrun
    def solutions(self, on_model, on_model_weight=None, limit=0,
                    force_weight=None):

        tsfile = self.termset.to_file()

        control = gringo.Control(["--conf=trendy", "--stats", "0", "--opt-strat=usc"])

        do_mincard = self.opts.family == "mincard" \
            or self.opts.force_size is not None
        do_subsets = self.opts.family == "subset" \
            or (self.opts.family =="mincard" and self.opts.mincard_tolerance)

        for f in self.domain:
            control.load(f)
        control.load(aspf("supportConsistency.lp"))
        control.load(aspf("minimizeWeightOnly.lp"))
        if do_mincard:
            control.load(aspf("minimizeSizeOnly.lp"))
        control.load(aspf("normalize.lp"))
        control.load(tsfile)

        control.ground([("base", [])])

        control.load(aspf("show.lp"))
        control.ground([("show", [])])

        start = time.time()

        if force_weight is None:
            control.assign_external(gringo.Fun("tolerance"),False)
            dbg("# start initial solving")
            opt = []
            res = control.solve(None, lambda model: opt.append(model.optimization()))
            dbg("# initial solve took %s" % (time.time()-start))

            optimizations = opt.pop()
            dbg("# optimizations = %s" % optimizations)

            weight = optimizations[0]
            if do_mincard:
                minsize = optimizations[1]
            if weight > 0 and on_model_weight is not None:
                for sample in self.solution_samples():
                    on_model_weight(sample)
                return

            control.assign_external(gringo.Fun("tolerance"),True)
        else:
            weight = force_weight
            dbg("# force weight = %d" % weight)

        max_weight = weight + self.opts.weight_tolerance
        control.add("minWeight", [], ":- not " + str(weight) + " #sum {Erg,E,T,S : measured(E,T,S,V), not guessed(E,T,S,V), toGuess(E,T,S), obs(E,T,S,M), Erg=50-M, M < 50;" + " Erg,E,T,S : measured(E,T,S,V), not guessed(E,T,S,V), toGuess(E,T,S), obs(E,T,S,M), Erg=M-49, M >= 50} " + str(max_weight) + " .")
        control.ground([("minWeight", [])])

        control.conf.solve.opt_mode = "ignore"
        control.conf.solve.project = 1 # ????
        control.conf.solve.models = limit # ????
        #print control.conf.solver[0].keys()

        if do_mincard:
            if self.opts.force_size:
                maxsize = self.opts.force_size
            else:
                maxsize = minsize + self.opts.mincard_tolerance
            control.add("minSize", [], ":- not " + str(minsize) + " #sum {L,I,J : dnf(I,J) , hyper(I,J,L)} " + str(maxsize) + ".")
            control.ground([("minSize", [])])

        if do_subsets:
            control.conf.solve.enum_mode = "domRec"
            control.conf.solver[0].heuristic = "Domain"
            control.conf.solver[0].dom_mod = "5,16"

        start = time.time()
        dbg("# begin enumeration")
        res = control.solve(None, on_model)
        dbg("# enumeration took %s" % (time.time()-start))


