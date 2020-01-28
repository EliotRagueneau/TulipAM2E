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
        self.initProperties()
        self.importInteractions()
        self.importGeneExpression()
        self.styleGraph()
        self.importPathways()
        return True

    def initProperties(self):
        self.idToNode = {}
        self.chromosome = self.graph.getStringProperty("chromosome")
        self.expression = self.graph.getStringProperty("expression")
        self.interactionStatus = self.graph.getStringProperty("interactionStatus")
        self.distance = self.graph.getDoubleProperty('distance')
        self.viewLayout = self.graph.getLayoutProperty('viewLayout')
        self.viewColor = self.graph.getColorProperty('viewColor')

    def importInteractions(self):
        dataFrame = read_csv(self.dataSet['Path to interaction csv'], sep='\t')
        for row in dataFrame.itertuples():
            id1, id2 = row.ID_locus1, row.ID_locus2
            for id in (id1, id2):
                if id not in self.idToNode:
                    self.idToNode[id] = self.graph.addNode({'viewLabel': id, 'chromosome': row.chromosome})
            self.graph.addEdge(self.idToNode[id1], self.idToNode[id2], {'distance': row.distance,
                                                                        'interactionStatus': row.interaction_status})

    def importGeneExpression(self):
        dataFrame = read_csv(self.dataSet['Path to expression csv'], sep='\t')
        for row in dataFrame.itertuples():
            try:
                self.expression[self.idToNode[row.IDs]] = str(row.expression)
            except KeyError:
                node = self.graph.addNode({'viewLabel': row.IDs, 'chromosome': row.chromosome})
                self.idToNode[row.IDs] = node
                self.expression[node] = str(row.expression)

    def importPathways(self):
        self.importPathwaysFromCSV(self.dataSet['Path to Reactome Symbols csv'])
        self.importPathwaysFromCSV(self.dataSet['Path to KEGG Symbols csv'])

    def importPathwaysFromCSV(self, csvPath):
        with open(csvPath, 'r') as csv:
            for line in csv:
                name, url, *loci = line.strip().split('\t')
                nodesOfPathway = []
                for locus in loci:
                    try:
                        nodesOfPathway.append(self.idToNode[locus])
                    except KeyError:
                        node = self.graph.addNode({'viewLabel': locus})
                        self.idToNode[locus] = node
                        nodesOfPathway.append(node)

                self.graph.inducedSubGraph(nodesOfPathway, self.graph, name)

    def styleGraph(self):
        self.applyLayout()
        self.colorNodes()
        self.colorEdges()

    def applyLayout(self):
        fm3Properties = tlp.getDefaultPluginParameters('FM^3 (OGDF)', graph=None)
        fm3Properties['Edge Length Property'] = self.distance
        self.graph.applyLayoutAlgorithm('FM^3 (OGDF)', self.viewLayout, fm3Properties)

    def colorNodes(self):
        pass

    def colorEdges(self):
        pass


# The line below does the magic to register the plugin into the plugin database
# and updates the GUI to make it accessible through the menus.
tulipplugins.registerPlugin("InteractionNetwork", "Interaction Network", "AM2E", "28/01/2020", "", "1.0")
