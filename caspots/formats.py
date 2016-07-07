
from zope import component, interface
from caspo import core

class ILogicalNetworkList(interface.Interface):
    """"""
    networks = interface.Attribute("")
    def __iter__(self):
        """"""
    def __len__(self):
        """"""

class LogicalNetworkList(list):
    interface.implements(ILogicalNetworkList)

    def __init__(self, networks=[], update_names=True):
        super(LogicalNetworkList, self).__init__(networks)
        if update_names:
            names = component.getUtility(core.ILogicalNames)
            for network in networks:
                names.add(network.mapping.itervalues())

    def add(self, network, update_names=True):
        super(LogicalNetworkList, self).add(network)
        if update_names:
            names = component.getUtility(core.ILogicalNames)
            names.add(network.mapping.itervalues())

class CsvReader2LogicalNetworkList(core.CsvReader2LogicalNetworkSet):
    interface.implements(ILogicalNetworkList)

    def __init__(self, reader):
        super(CsvReader2LogicalNetworkList, self).__init__(reader)

    def _read(self, reader):
        self._networks = LogicalNetworkList(map(lambda row: core.ILogicalNetwork(core.LogicalMapping(row)), reader))

gsm = component.getGlobalSiteManager()
gsm.registerAdapter(CsvReader2LogicalNetworkList, (core.ICsvReader,), ILogicalNetworkList)

