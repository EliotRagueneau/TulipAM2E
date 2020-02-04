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
        self.pluginProgress.progress(0, 4)
        self.initProperties()
        self.importInteractions()
        self.pluginProgress.progress(1, 4)
        self.importGeneExpression()
        self.pluginProgress.progress(2, 4)
        self.importPathways()
        self.pluginProgress.progress(3, 4)
        self.styleGraph()
        return True

    def initProperties(self):
        self.idToNode = {}
        self.chromosome = self.graph.getStringProperty("chromosome")
        self.expression = self.graph.getStringProperty("expression")
        self.interactionStatus = self.graph.getStringProperty("interactionStatus")
        self.distance = self.graph.getIntegerProperty('distance')
        self.viewLayout = self.graph.getLayoutProperty('viewLayout')
        self.viewColor = self.graph.getColorProperty('viewColor')
        self.viewSize = self.graph.getSizeProperty('viewSize')

    def importInteractions(self):
        dataFrame = read_csv(self.dataSet['Path to interaction csv'], sep='\t')
        for row in dataFrame.itertuples():
            chromosome = row.chromosome
            id1, id2 = row.ID_locus1, row.ID_locus2
            for id in (id1, id2):
                if id not in self.idToNode:
                    label = "" if "chr6" in id else id
                    self.idToNode[id] = self.graph.addNode({'viewLabel': label, 'chromosome': chromosome})
            self.graph.addEdge(self.idToNode[id1], self.idToNode[id2], {'distance': int(row.distance),
                                                                        'interactionStatus': row.interaction_status,
                                                                        'chromosome': chromosome})

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
                        gene = self.idToNode[locus]
                        nodesOfPathway.append(gene)
                        for neigh in self.graph.getInOutNodes(gene):
                            nodesOfPathway.append(neigh)
                    except KeyError:
                        node = self.graph.addNode({'viewLabel': locus})
                        self.idToNode[locus] = node
                        nodesOfPathway.append(node)

                self.graph.inducedSubGraph(nodesOfPathway, self.graph, name)

    def styleGraph(self):
        self.applyLayout()
        self.colorLoci()
        self.colorInteractions()
        # self.setLociSize()

    def applyLayout(self):
        fm3Properties = tlp.getDefaultPluginParameters('FM^3 (OGDF)', graph=None)
        fm3Properties['Edge Length Property'] = self.distance
        self.graph.applyLayoutAlgorithm('FM^3 (OGDF)', self.viewLayout, fm3Properties)

    def colorLoci(self):
        self.colorLociByExpression('intergenic', tlp.Color(222, 212, 195, 80))
        self.colorLociByExpression('up', tlp.Color(0, 255, 0))
        self.colorLociByExpression('down', tlp.Color(255, 0, 0))
        self.colorLociByExpression('stable', tlp.Color(0, 0, 0, 80))
        self.colorLociByExpression('', tlp.Color(0, 0, 0, 80))
        self.colorLociByExpression('nan', tlp.Color(255, 230, 255, 80))

    def colorLociByExpression(self, expression, color):
        for loci in self.expression.getNodesEqualTo(expression):
            self.viewColor[loci] = color

    def colorInteractions(self):
        self.colorInteractionsByStatus('gain', tlp.Color(0, 255, 0))
        self.colorInteractionsByStatus('stable', tlp.Color(0, 0, 0, 80))
        self.colorInteractionsByStatus('loss', tlp.Color(255, 0, 0))

    def colorInteractionsByStatus(self, interactionStatus, color):
        for interaction in self.interactionStatus.getEdgesEqualTo(interactionStatus):
            self.viewColor[interaction] = color

    def setLociSize(self):
        self.viewSize.setAllNodeValue(tlp.Vec3f(1e5, 1e5, 0.5))


tulipplugins.registerPlugin("InteractionNetwork", "Interaction Network", "AM2E", "28/01/2020", "", "1.0")
