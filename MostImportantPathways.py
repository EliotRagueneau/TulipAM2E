import tulipplugins
from tulip import tlp


class ImportantPathways(tlp.Algorithm):
    def __init__(self, context):
        tlp.Algorithm.__init__(self, context)

    def initProperties(self):
        self.expression = self.graph.getStringProperty("expression")
        self.chromosome = self.graph.getStringProperty("chromosome")
        self.interactionStatus = self.graph.getStringProperty("interactionStatus")

    def run(self):
        self.initProperties()
        pathwaysRatios = []
        modalities = ('up', 'down', 'gain', 'loss', 'ch6')
        for pathway in self.graph.getSubGraphs():
            ratios = {modality: 0 for modality in modalities}
            self.evaluatePathwayLociRatios(pathway, ratios)
            self.evaluatePathwayInteractionsRatios(pathway, ratios)
            pathwaysRatios.append((pathway.getName(), ratios))

        self.filterAndSortRatiosByModality(pathwaysRatios, modalities)
        return True

    def evaluatePathwayLociRatios(self, pathway, ratios):
        n = pathway.numberOfNodes()

        for locus in pathway.getNodes():
            for expressionValue in ('up', 'down'):
                if self.expression[locus] == expressionValue:
                    ratios[expressionValue] += 1
            if self.chromosome[locus] == 'chr6':
                ratios['ch6'] += 1

        for param in ('up', 'down'):
            if ratios['ch6'] == 0:
                ratios[param] = 0
            else:
                ratios[param] /= ratios['ch6']
        ratios['ch6'] /= n

    def evaluatePathwayInteractionsRatios(self, pathway, ratios):
        e = pathway.numberOfEdges()
        for interactionStatus in ('gain', 'loss'):
            if e == 0:
                ratios[interactionStatus] = 0
            else:
                for interaction in pathway.getEdges():
                    if self.interactionStatus[interaction] == interactionStatus:
                        ratios[interactionStatus] += 1
                ratios[interactionStatus] /= e

    def filterAndSortRatiosByModality(self, pathwaysRatios, modalities):
        for modality in modalities:
            pathwaysRatios.sort(key=lambda elt: (elt[1][modality], elt[1]['ch6']), reverse=True)
            print('\n', modality.capitalize())
            print('name,', ','.join(modalities))
            for elt in pathwaysRatios:
                if elt[1][modality] != 0:
                    print('{}, '.format(elt[0]), end='')
                    for param in modalities:
                        print('{:.2%}, '.format(elt[1][param]), end='')
                    print()


tulipplugins.registerPluginOfGroup("ImportantPathways", "Important Pathways", "AM2E",
                                   "29/01/2020", "", "1.0", "AM2E")
