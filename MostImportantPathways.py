# When the plugin development is finished, you can copy the associated
# Python file to /net/cremi/edarnige/.Tulip-5.4/plugins/python
# or /autofs/unityaccount/ens/tulip/lib/tulip/python/
# and it will be automatically loaded at Tulip startup
import tulipplugins
from tulip import tlp


class ImportantPathways(tlp.Algorithm):
    def __init__(self, context):
        tlp.Algorithm.__init__(self, context)
        # You can add parameters to the plugin here through the
        # following syntax:
        # self.add<Type>Parameter('<paramName>', '<paramDoc>',
        #                         '<paramDefaultValue>')
        # (see the documentation of class tlp.WithParameter to see what
        #  parameter types are supported).

    def initProperties(self):
        self.expression = self.graph.getStringProperty("expression")
        self.chromosome = self.graph.getStringProperty("chromosome")
        self.interactionStatus = self.graph.getStringProperty(
            "interactionStatus")

    def check(self):
        # This method is called before applying the algorithm on the
        # input graph. You can perform some precondition checks here.
        # See comments in the run method to know how to have access to
        # the input graph.
        #
        # Must return a tuple (Boolean, string). First member indicates if the
        # algorithm can be applied and the second one can be used to provide
        # an error message.
        return (True, '')

    def run(self):
        # This method is the entry point of the algorithm when it is called
        # and must contain its implementation.
        #
        # The graph on which the algorithm is applied can be accessed through
        # the 'graph' class attribute (see documentation of class tlp.Graph).
        #
        # The parameters provided by the user are stored in a dictionary
        # that can be accessed through the 'dataSet' class attribute.
        #
        # The method must return a Boolean indicating if the algorithm
        # has been successfully applied on the input graph.
        self.initProperties()
        data = []
        modalities = ('up', 'down', 'gain', 'loss', 'ch6')
        for sub in self.graph.getSubGraphs():
            params = {modality: 0 for modality in modalities}
            n = sub.numberOfNodes()
            e = sub.numberOfEdges()

            for locus in sub.getNodes():
                for expressionValue in ('up', 'down'):
                    if self.expression[locus] == expressionValue:
                        params[expressionValue] += 1
                if self.chromosome[locus] == 'chr6':
                    params['ch6'] += 1

            for param in ('up', 'down'):
                if params['ch6'] == 0:
                    params[param] = 0
                else:
                    params[param] /= params['ch6']
            params['ch6'] /= n

            for interactionStatus in ('gain', 'loss'):
                if e == 0:
                    params[interactionStatus] = 0
                else:
                    for interaction in sub.getEdges():
                        if self.interactionStatus[interaction] == interactionStatus:
                            params[interactionStatus] += 1
                    params[interactionStatus] /= e
            data.append((sub.getName(), params))

        for modality in modalities:
            data.sort(key=lambda elt: (elt[1][modality], elt[1]['ch6']), reverse=True)
            print('\n', modality.capitalize())
            print('name,', ','.join(modalities))
            for elt in data:
                if elt[1][modality] != 0:
                    print('{}, '.format(elt[0]), end='')
                    for param in modalities:
                        print('{:.2%}, '.format(elt[1][param]), end='')
                    print()

        return True


# The line below does the magic to register the plugin into the plugin database
# and updates the GUI to make it accessible through the menus.
tulipplugins.registerPluginOfGroup(pluginClassName='ImportantPathways',
                                   pluginName='Important Pathways',
                                   author='AM2E',
                                   date='29/01/2020',
                                   info='find most important pathways based on up/down regulating loci',
                                   release='1.0',
                                   group='AM2E')
