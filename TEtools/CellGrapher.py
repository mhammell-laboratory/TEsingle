import networkx as nx


def wildcard_adder(wdict, umi):
    match_list = []
    for i in range(len(umi)):
        wildumi = umi[:i]+'*'+umi[i+1:]
        if wildumi in wdict:
            match_list += wdict[wildumi]
            wdict[wildumi].append(umi)
        else:
            wdict[wildumi] = [umi]

    return wdict, set(match_list)


class cell_graph(object):
    def __init__(self):
        self.G = nx.Graph()
        self.subgraphs = []

    def grapher(self, umidict):
        wilddict = {}
        for umi in umidict:
            wilddict, match_set = wildcard_adder(wilddict, umi)
            self.G.add_node(umi)
            for mumi in match_set:
                self.G.add_edge(umi, mumi)

        del wilddict
        return self.G

    def subgrapher(self):
        self.subgraphs = [self.G.subgraph(c).copy() for c in nx.connected_components(self.G)]
        return self.subgraphs

    def build(self, umidict):
        self.G = self.grapher(umidict)
        self.subgraphs = self.subgrapher()