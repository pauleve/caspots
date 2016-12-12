
import sys

import pandas as pd
import numpy as np

import gringo

from caspo.core.setup import Setup
from caspo.core.literal import Literal
from caspo.core.clamping import Clamping, ClampingList

from asputils import *


class Experiment:
    def __init__(self, id):
        self.id = id
        self.obs = {}
        self.dobs = {}
        self.mutations = {}

    def add_mutation(self, node, clamp):
        self.mutations[node] = clamp

    def add_obs(self, t, node, value, dvalue):
        if t not in self.obs:
            self.obs[t] = {}
            self.dobs[t] = {}
        self.obs[t][node] = value
        self.dobs[t][node] = dvalue

    def __str__(self):
        buf = "Experiment(%d):\n" % self.id
        for a, c in sorted(self.mutations.items()):
            buf += "\t%2d %s\n" % (c,a)
        buf += "----\n"
        for t, values in sorted(self.obs.items()):
            buf += "\t%4d |" % t
            for n, v in sorted(values.items()):
                buf += "\t%s=%d" % (n,v)
            buf += "\n"
        buf += "----"
        return buf

class Dataset:
    def __init__(self, name, dfactor=100, discretize="round"):
        self.name = name
        self.stimulus = set()
        self.inhibitors = set()
        self.readout = set()
        self.experiments = {}
        self.dfactor = dfactor
        self.discretize = getattr(self, "discretize_%s" % discretize)

    def discretize_round(self, value):
        return int(round(self.dfactor*value))

    def binarize(self, dvalue):
        return 1 if dvalue >= self.dfactor/2 else 0

    def load_from_midas(self, midas, graph):
        df = pd.read_csv(midas)
        df.drop(df.columns[0], axis=1, inplace=True)
        df = df.reset_index(drop=True)

        def is_stimulus(name):
            return name.startswith('TR') and not name.endswith('i')
        def is_inhibitor(name):
            return name.startswith('TR') and name.endswith('i')
        def is_readout(name):
            return name.startswith('DV')

        stimuli = [c[3:] for c in [c for c in df.columns if is_stimulus(c)]]
        inhibitors = [c[3:-1] for c in [c for c in df.columns if is_inhibitor(c)]]
        readouts = [c[3:] for c in [c for c in df.columns if is_readout(c)]]
        self.setup = Setup(stimuli, inhibitors, readouts)

        exp_t = {}
        def exp_of_clamps(clamps):
            if clamps not in exp_t:
                eid = len(exp_t)
                exp = Experiment(eid)
                for node, clamp in clamps:
                    exp.add_mutation(node, clamp)
                self.experiments[eid] = exp
                exp_t[clamps] = exp
                return exp
            else:
                return exp_t[clamps]

        for i, row in df.iterrows():
            clamps = set()
            for var, sign in row.filter(regex='^TR').iteritems():
                var = var[3:]
                sign = int(sign)
                if var in stimuli:
                    if sign == 1 or len(graph.predecessors(var)):
                        clamps.add((var, sign or -1))
                elif sign == 1:
                    clamps.add((var, -1))
            clamps = tuple(clamps)
            exp = exp_of_clamps(clamps)

            for var, fvalue in row.filter(regex='^DV').iteritems():
                if np.isnan(fvalue):
                    continue
                var = var[3:]
                time = int(row.get("DA:%s" % var))
                dvalue = self.discretize(fvalue)
                bvalue = self.binarize(dvalue)
                exp.add_obs(time, var, bvalue, dvalue)

    def to_funset(self):
        fs = funset(self.setup)

        clampings = []
        for exp in sorted(self.experiments.values(), key=lambda e: e.id):
            i = exp.id
            literals = [Literal(node, sign) for node, sign in \
                            exp.mutations.items()]
            clampings.append(Clamping(literals))
            for time, obs in exp.dobs.items():
                for var, dval in obs.items():
                    fs.add(gringo.Fun('obs', [i, time, var, dval]))
        clampings = ClampingList(clampings)
        fs.update(clampings.to_funset("exp"))

        fs.add(gringo.Fun('dfactor', [self.dfactor]))
        return fs

    """
    def feed_from_asp(self, fs):
        def key(f):
            if f.name() == "exp":
                return 0
            return 10
        for f in sorted(fs, key=key):
            p = f.name()
            args = f.args()
            if p == "stimulus":
                self.stimulus.add(args[0])
            elif p == "inhibitor":
                self.inhibitors.add(args[0])
            elif p == "readout":
                self.readout.add(args[0])
            elif p == "exp":
                id = args[0]
                self.experiments[id] = Experiment(id)
            elif p == "dfactor":
                assert args[0] == 100, args
                continue
            elif p == "clamped":
                eid, node, clamp = args
                self.experiments[eid].add_mutation(node, clamp)
            elif p == "obs":
                eid, t, node, value = args
                value = 1 if value >= 50 else 0
                self.experiments[eid].add_obs(t, node, value)
            else:
                raise Exception("unknown clause '%s'" % p)
    """

    def __str__(self):
        buf = "%s %s %s\n" % ("#"*10, self.name, "#"*10)
        buf += "\n".join(map(str, self.experiments.values()))
        return buf

