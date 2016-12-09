

import caspo.learn
import caspo.core

from caspots import asputils
from .model import *

"""
class TimeSeries2TermSet(asp.TermSetAdapter):
    component.adapts(core.IDataset, learn.IDiscretization)
    def __init__(self, dataset, discretization):
        super(TimeSeries2TermSet, self).__init__()
        self._termset = self._termset.union(asp.interfaces.ITermSet(dataset.setup))
        self._termset.add(asp.Term('dfactor', [discretization.factor]))
        for i, cues in enumerate(dataset.cues):
            self._termset.add(asp.Term('exp', [i]))
            self._termset = self._termset.union(component.getMultiAdapter((cues, dataset), asp.ITermSet))
            for time in dataset.times:
                for name, value in dataset.obs[i][time].iteritems():
                    self._termset.add(asp.Term('obs', [i, time, name, discretization(value)]))
"""

class Csv2Dataset(caspo.core.Dataset):
    def __init__(self, midas, graph):
        df = pd.read_csv(midas)
        self.graph = graph
        self.times = np.unique(df.filter(regex='^DA').values.flatten())
        df.drop(df.columns[0], axis=1, inplace=True)
        super(caspo.Dataset, self).__init__(df.reset_index(drop=True))

        stimuli = [c[3:] for c in [c for c in self.columns if self.is_stimulus(c)]]
        inhibitors = [c[3:-1] for c in [c for c in self.columns if self.is_inhibitor(c)]]
        readouts = [c[3:] for c in [c for c in self.columns if self.is_readout(c)]]

        self.setup = Setup(stimuli, inhibitors, readouts)

    @property
    def clampings(self):
        clampings = []
        for _, row in self.filter(regex='^TR').iterrows():
            literals = []
            for var, sign in row.iteritems():
                if self.is_stimulus(var):
                    node = var[3:]
                    # this is different from caspo: negative literal only if
                    # parents
                    if sign == 1:
                        literals.append(Literal(node, 1))
                    elif len(self.graph.predecessors(node)):
                        literals.append(Literal(node, -1))
                else:
                    if sign == 1:
                        literals.append(Literal(var[3:-1], -1))

            clampings.append(Clamping(literals))

        return ClampingList(clampings)

    """
    def at(self, time):
        if time not in self.times:
            raise ValueError("The time-point %s does not exists in the dataset. Available time-points are: %s" % (time, list(self.times)))

        cond = True
        for col in df.filter(regex='^DA').columns:
            cond = cond & (self[col] == time)

        for i, (cues, obs) in enumerate(zip(self.cues, self.obs)):
            yield i, cues, obs[time]
    """

def parse_answer(line):
    return asputils.re_answer.findall(line)

def termset_of_sif(siffile):
    sif = component.getUtility(core.IFileReader)
    sif.read(siffile)

    graph = core.IGraph(sif)
    length = core.Length(0)

    names = component.getUtility(core.ILogicalNames)
    names.load(graph, length.length)
    termset = asp.ITermSet(names)
    return termset, graph.graph

def termset_of_dataset(csvfile, graph, discretization="round", factor=100):
    dataset = Csv2Dataset(csvfile, graph)
    discretize = partial(getattr(caspo.learn.Learner, discretization), factor)
    fs = funset()
    fs.update(dataset.to_funset(discretize))
    fs.add(gringo.Fun('dfactor', [factor]))
    return fs

def domain_of_networks(networks, hypergraph, dataset):
    out = ["1{%s}1." % ("; ".join(["model(%d)" % i for i in range(len(networks))]))]

    nodefid = {}
    edgeids = {}
    hids = {}
    for t in funset(hypergraph):
        p = t.pred
        if p == "node":
            nodefid[t.arg(0)] = t.arg(1)
        elif p == "edge":
            hid, node, sign = (t.arg(0), t.arg(1), t.arg(2))
            k = (node, sign)
            if k not in edgeids:
                edgeids[k] = set()
            edgeids[k].add(hid)
        elif p == "hyper":
            fid, eid, n = (t.arg(0), t.arg(1), t.arg(2))
            k = (fid,n)
            if k not in hids:
                hids[k] = set()
            hids[k].add(eid)

    nodes = dataset.inhibitors.union(dataset.readout)

    predicates = ["formula", "dnf", "clause"]

    # {formula(V,I): node(V,I)}6, 6{dnf(I,J): hyper(I,J,N)}6, 7{clause(J,V,B): edge(J,V,B)}7

    #convert and write each network to logic facts
    for i, net in enumerate(networks):
        ts = asp.ITermSet(net)
        data = ts.to_str().replace(".\n", " ")
        data = parse_answer(data)

        bn = BNModel(0)
        bn.feed_from_asp(data)
        bn.default_false(nodes)

        facts = []
        for f in bn.formula.values():
            fid = nodefid[f.var]
            facts.append("formula(\"%s\",%d)" % (f.var, fid))
            for dnf in f:
                hid = None
                n = len(dnf)
                for (node, sign) in dnf:
                    eids = edgeids[(node,sign)].copy()
                    if hid is None:
                        hid = eids
                    else:
                        hid.intersection_update(eids)
                hid.intersection_update(hids[(fid,n)])
                assert len(hid) == 1, hid
                hid = hid.pop()
                facts.append("dnf(%d,%d)" % (fid, hid))
                for (node, sign) in dnf:
                    facts.append("clause(%d,\"%s\",%d)" % (hid, node, sign))

        for f in facts:
            out.append("%s :- model(%d)." % (f,i))
    return "\n".join(out) + "\n"


