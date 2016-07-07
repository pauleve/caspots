
import sys

from asputils import re_clause, parse_args


class Experiment:
    def __init__(self, id):
        self.id = id
        self.obs = {}
        self.mutations = {}

    def add_mutation(self, node, clamp):
        self.mutations[node] = clamp

    def add_obs(self, t, node, value):
        if t not in self.obs:
            self.obs[t] = {}
        self.obs[t][node] = value

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
    def __init__(self, name):
        self.name = name
        self.stimulus = set()
        self.inhibitors = set()
        self.readout = set()
        self.experiments = {}

    def feed_from_asp(self, content):
        def key(e):
            try:
                return ["exp"].index(e[0])
            except ValueError:
                return 10
        for (p, args) in sorted(re_clause.findall(content), key=key) :
            try:
                args = parse_args(args)
            except ValueError:
                continue
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

    def __str__(self):
        buf = "%s %s %s\n" % ("#"*10, self.name, "#"*10)
        buf += "\n".join(map(str, self.experiments.values()))
        return buf


if __name__ == "__main__":
    import os
    import sys
    lpfile = sys.argv[1]
    name = os.path.basename(lpfile).replace(".lp", "")
    d = Dataset(name)
    d.feed_from_asp(open(lpfile).read())
    print(str(d))


