
import math
import sys
from collections import defaultdict
from zope import component, interface
from pyzcasp import asp
from caspo import core, learn
from caspo.core import impl

from caspots import asputils
from caspots.model import *

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



class CsvReader2Dataset(object):
    component.adapts(core.ICsvReader, core.IGraph)
    interface.implements(core.IDataset, core.IClampingList)

    def __init__(self, reader, graph):
        #Header without the CellLine column
        species = reader.fieldnames[1:]    
        stimuli = map(lambda name: name[3:], filter(lambda name: name.startswith('TR:') and not name.endswith('i'), species))
        inhibitors = map(lambda name: name[3:-1], filter(lambda name: name.startswith('TR:') and name.endswith('i'), species))
        readouts = map(lambda name: name[3:], filter(lambda name: name.startswith('DV:'), species))

        self.setup = impl.Setup(stimuli, inhibitors, readouts)

        self.cues = []
        self.obs = []
        self.nobs = defaultdict(int)
                
        times = []
        for row in reader:
            literals = []
            for s in stimuli:
                if row['TR:' + s] == '1':
                    literals.append(impl.Literal(s,1))
                elif not list(graph.predecessors(s)):
                    literals.append(impl.Literal(s,-1))
                    
            for i in inhibitors:
                if row['TR:' + i + 'i'] == '1':
                    literals.append(impl.Literal(i,-1))

            clamping = impl.Clamping(literals)
            obs = defaultdict(dict)
            for r in readouts:
                if not math.isnan(float(row['DV:' + r])):
                    time = int(row['DA:' + r]) 
                    times.append(time)
                    obs[time][r] = float(row['DV:' + r]) 
                    self.nobs[time] += 1
                    
            if clamping in self.cues:
                index = self.cues.index(clamping)
                self.obs[index].update(obs)
            else:
                self.cues.append(clamping)
                self.obs.append(obs)
                
        self.times = frozenset(times)
        self.nexps = len(self.cues)
        
    @property
    def clampings(self):
        return self.cues
        
    def at(self, time):
        if time not in self.times:
            raise ValueError("The time-point %s does not exists in the dataset. Available time-points are: %s" % (time, list(self.times)))
                                  
        for i, (cues, obs) in enumerate(zip(self.cues, self.obs)):
            yield i, cues, obs[time]


gsm = component.getGlobalSiteManager()
gsm.registerAdapter(TimeSeries2TermSet)
gsm.registerAdapter(CsvReader2Dataset, (core.ICsvReader,core.IGraph), core.IDataset)

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
    reader = component.getUtility(core.ICsvReader)
    reader.read(csvfile)
    dataset = component.getMultiAdapter((reader, graph), core.IDataset)
    #dataset = core.IDataset(reader)
    discretize = component.createObject(discretization, factor)
    return component.getMultiAdapter((dataset, discretize), asp.ITermSet)

def domain_of_networks(networks, pkn, dataset):
    out = ["1{%s}1." % ("; ".join(["model(%d)" % i for i in range(len(networks))]))]

    nodefid = {}
    edgeids = {}
    hids = {}
    for t in pkn:
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


