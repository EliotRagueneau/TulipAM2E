# When the plugin development is finished, you can copy the associated Python file 
# to /home/eliot/.Tulip-5.3/plugins/python
# or /usr/local/lib/tulip/python/
# and it will be automatically loaded at Tulip startup

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

    def importGraph(self):
        self.import_interactions()
        self.import_gene_expression()

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
            if line.IDs not in self.nodes:
                self.nodes[line.IDs] = self.graph.addNode()
            expression[self.nodes[line.IDs]] = str(line.expression)

# The line below does the magic to register the plugin into the plugin database
# and updates the GUI to make it accessible through the menus.
tulipplugins.registerPlugin("InteractionNetwork", "Interaction Network", "AM2E", "28/01/2020", "", "1.0")
