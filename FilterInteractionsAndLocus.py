import tulipplugins
from tulip import tlp


class FilterInteractionsAndLocus(tlp.Algorithm):
    def __init__(self, context):
        tlp.Algorithm.__init__(self, context)

    def notify(self, step, max_step, message):
        self.pluginProgress.setComment(message)
        self.pluginProgress.progress(step, max_step)

    def run(self):
        self.notify(0, 4, "Initialize parameters")
        self.initProperties()
        self.notify(1, 4, "Filter Loci and Interactions by their expression and status")
        self.filterLociAndInteractions()
        self.notify(2, 4, "Filter empty pathways")
        self.filterPathways()
        self.notify(3, 4, "Apply layout algorithm")
        self.styleGraph()
        self.pluginProgress.progress(4, 4)
        return True

    def initProperties(self):
        self.expression = self.graph.getStringProperty('expression')
        self.interactionStatus = self.graph.getStringProperty('interactionStatus')
        self.chromosome = self.graph.getStringProperty("chromosome")
        self.viewLayout = self.graph.getLocalLayoutProperty('viewLayout')
        self.viewSize = self.graph.getSizeProperty('viewSize')

    def filterLociAndInteractions(self):
        toDelete = [locus for locus in self.graph.getNodes() if not self.isLocusRelevant(locus, self.graph, False)]
        self.graph.delNodes(toDelete, True)

    def isLocusRelevant(self, locus, graph, is_neigh=False):
        if self.chromosome[locus] == '' or self.chromosome[locus] is None:
            return False

        if self.expression[locus] in ('up', 'down'):
            return True

        for neigh in graph.getInOutNodes(locus):
            if not is_neigh and self.isLocusRelevant(neigh, graph, True):
                return True

        for interaction in graph.getInOutEdges(locus):
            if self.isInteractionRelevant(interaction):
                return True

        return False

    def isInteractionRelevant(self, interaction):
        return self.interactionStatus[interaction] in ('gain', 'loss')

    def filterPathways(self):
        toDelete = [pathway for pathway in self.graph.getSubGraphs() if not self.isPathwayRelevant(pathway)]
        for pathway in toDelete:
            self.graph.delSubGraph(pathway)

    def isPathwayRelevant(self, pathway):
        for locus in pathway.getNodes():
            if self.graph.isElement(locus):
                if self.chromosome[locus] and self.expression[locus] in ('up', 'down'):
                    return True
        return False

    def styleGraph(self):
        self.applyLayout()
        self.setLociSize()

    def applyLayout(self):
        fm3Properties = tlp.getDefaultPluginParameters('FM^3 (OGDF)', graph=None)
        self.graph.applyLayoutAlgorithm('FM^3 (OGDF)', self.viewLayout, fm3Properties)

    def setLociSize(self):
        self.viewSize.setAllNodeValue(tlp.Vec3f(1e4, 1e4, 0.5))

tulipplugins.registerPluginOfGroup("FilterInteractionsAndLocus", "Filter interactions and locus", "AM2E",
                                   "29/01/2020", "", "1.0", "AM2E")
