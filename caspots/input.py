
import math
import sys
from collections import defaultdict
from zope import component, interface
from pyzcasp import asp
from caspo import core, learn
from caspo.core import impl

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

