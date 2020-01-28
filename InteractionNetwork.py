import tulipplugins
from pandas import read_csv
from tulip import tlp


class InteractionNetwork(tlp.ImportModule):
    def __init__(self, context):
        tlp.ImportModule.__init__(self, context)
        self.addFileParameter("Path to interaction csv", help="",
                              defaultValue="/home/eliot/Documents/Travail/M2/DEA/Tulip/ProjetM2/interactions_chromosome6.csv",
                              isMandatory=True, mustExist=True)
        self.addFileParameter("Path to expression csv", help="",
                              defaultValue="/home/eliot/Documents/Travail/M2/DEA/Tulip/ProjetM2/chromosome6_fragments_expressions.csv",
                              isMandatory=True, mustExist=True)
        self.addFileParameter("Path to Reactome Symbols csv", help="",
                              defaultValue="/home/eliot/Documents/Travail/M2/DEA/Tulip/ProjetM2/REACTOME.symbols.csv",
                              isMandatory=True, mustExist=True)
        self.addFileParameter("Path to KEGG Symbols csv", help="",
                              defaultValue="/home/eliot/Documents/Travail/M2/DEA/Tulip/ProjetM2/KEGG.symbols.csv",
                              isMandatory=True, mustExist=True)

    def importGraph(self):
        self.import_interactions()
        self.import_gene_expression()
        self.style_graph()
        self.import_pathways()

        return True

    def import_interactions(self):
        data_frame = read_csv(self.dataSet['Path to interaction csv'], sep='\t')
        self.nodes = {}
        for line in data_frame.itertuples():
            id1, id2 = line.ID_locus1, line.ID_locus2
            for id in (id1, id2):
                if id not in self.nodes:
                    self.nodes[id] = self.graph.addNode()
            self.graph.addEdge(self.nodes[id1], self.nodes[id2], {'distance': line.distance,
                                                                  'interaction_status': line.interaction_status})

    def import_gene_expression(self):
        data_frame = read_csv(self.dataSet['Path to expression csv'], sep='\t')
        expression = self.graph.getStringProperty("expression")
        for line in data_frame.itertuples():
            try:
                expression[self.nodes[line.IDs]] = str(line.expression)
            except KeyError:
                node = self.graph.addNode()
                self.nodes[line.IDs] = node
                expression[node] = str(line.expression)

    def import_pathways(self):
        self.import_pathways_from_csv(self.dataSet['Path to Reactome Symbols csv'])
        self.import_pathways_from_csv(self.dataSet['Path to KEGG Symbols csv'])

    def import_pathways_from_csv(self, csv_path):
        with open(csv_path, 'r') as csv:
            for line in csv:
                name, url, *loci = line.strip().split('\t')
                nodes_of_pathway = []
                for locus in loci:
                    try:
                        nodes_of_pathway.append(self.nodes[locus])
                    except KeyError:
                        node = self.graph.addNode()
                        self.nodes[locus] = node
                        nodes_of_pathway.append(node)

                self.graph.inducedSubGraph(nodes_of_pathway, self.graph, name)

    def style_graph(self):
        pass


# The line below does the magic to register the plugin into the plugin database
# and updates the GUI to make it accessible through the menus.
tulipplugins.registerPlugin("InteractionNetwork", "Interaction Network", "AM2E", "28/01/2020", "", "1.0")
