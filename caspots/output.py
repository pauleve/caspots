
import sys
import os
from zope import component
from pyzcasp import asp
from caspo import core, learn

def csv_of_resultlp(lpfile, pknfile, csvfile):

    sif = component.getUtility(core.IFileReader)
    sif.read(pknfile)
    graph = core.IGraph(sif)

    names = component.getUtility(core.ILogicalNames)
    names.load(graph)

    models = []
    next_is_model = False
    with open(lpfile) as f:
        for line in f:
            if next_is_model:
                next_is_model = False
                models.append(line)
            elif line.startswith("Answer"):
                next_is_model = True

    def network_of_asp(model):
        answer = asp.AnswerSet(model.split())
        return core.ILogicalNetwork(asp.ITermSet(answer))

    networks = list(map(network_of_asp, models))
    networks = core.LogicalNetworkSet(networks)

    writer = core.ICsvWriter(networks)
    writer.write(os.path.basename(csvfile), os.path.dirname(csvfile))

