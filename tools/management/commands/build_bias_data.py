from decimal import Decimal
import logging
import math

from django.db import models
from build.management.commands.base_build import Command as BaseBuild
from protein.models import ProteinGProteinPair
from ligand.models import *
from common.models import Publication




MISSING_PROTEINS = {}
SKIPPED = 0

class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)

    help = 'Reads bias data and imports it'
    # source file directory
    # structure_data_dir = os.sep.join([settings.EXCEL_DATA, 'ligand_data', 'bias'])
    structure_data_dir = 'excel/'
    publication_cache = {}
    ligand_cache = {}
    data_all = []

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-f', '--filename',
                            action='append',
                            dest='filename',
                            help='Filename to import. Can be used multiple times')
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing bias records')
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run',
                            default=False)

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        # delete any existing structure data
        if options['purge']:
            try:
                print('Started purging bias data')
                self.purge_bias_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
        # import the structure data
        try:
            print('CREATING BIAS DATA')
            print(options['filename'])
            self.build_bias_data()
            self.build_bias_data_subtypes()
            self.logger.info('COMPLETED CREATING BIAS')
        except Exception as msg:
            print('--error--', msg, '\n')

# pylint: disable=R0201
    def purge_bias_data(self):
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()

    # pylint: disable=C0321
    def get_from_model(self):
        try:
            content = BiasedExperiment.objects.all().prefetch_related(
                'experiment_data', 'ligand', 'receptor', 'publication'
                , 'publication__web_link'
                , 'experiment_data__emax_ligand_reference',
            ).order_by('publication', 'receptor', 'ligand')
        except:
            print('no data returned')
            content = None
        return content

    def process_data(self, content):
    # Disable all the no-member violations in this function
    # pylint: disable=no-member
        '''
        Merge BiasedExperiment with its children
        and pass it back to loop through dublicates
        '''
        rd = []
        counter = 0
        for instance in enumerate(content):
            temp_obj = []
            fin_obj = {}
            fin_obj['main'] = (instance[1])
            vendor_counter = 0
            for i in instance[1].experiment_data_vendors.all():
                vendor_counter = vendor_counter + 1
            for entry in instance[1].experiment_data.all():
                author_list = list()
                for author in entry.experiment_data_authors.all():
                    author_list.append(author.author)
                temp_obj.append(entry)
                counter += 1
            fin_obj['authors'] = author_list
            fin_obj['children'] = temp_obj
            fin_obj['vendor_counter'] = vendor_counter
            rd.append(fin_obj)
        return rd

    def queryset_to_dict(self, results):
        '''
        Merge bias experminet data with assay data
        '''
        send = list()
        for j in results:
            temp_dict = dict()
            temp = dict()
            doubles = []
            temp['publication'] = j['main'].publication
            temp['species'] = j['main'].receptor.species.common_name
            # temp['ligand'] = j['main'].ligand
            temp['endogenous_ligand'] = j['main'].endogenous_ligand
            temp['receptor'] = j['main'].receptor
            temp['assay'] = dict()
            temp['vendor_counter'] = j['vendor_counter']
            temp['reference'] = list()
            temp['authors'] = j['authors']
            temp['ref_ligand_experiment'] = dict()
            temp['article_quantity'] = 0
            temp['labs_quantity'] = 0
            temp['reference_ligand'] = None
            if not j['children']:
                continue
            temp_dict['potency'] = ''
            temp_dict['t_factor'] = ''
            temp_dict['log_bias_factor'] = ''
            temp_dict['order_no'] = 0
            temp_dict['reference_ligand'] = None
            temp_dict['signalling_protein'] = j['children'][0].signalling_protein.lower()
            temp_dict['cell_line'] = j['children'][0].cell_line
            temp_dict['family'] = j['children'][0].family
            temp_dict['assay_type'] = j['children'][0].assay_type
            temp_dict['assay_measure_method'] = j['children'][0].assay_measure
            temp_dict['assay_time_resolved'] = j['children'][0].assay_time_resolved
            if j['children'][0].quantitive_activity:
                temp_dict['quantitive_activity'] = j['children'][0].quantitive_activity
                temp_dict['quantitive_activity_initial'] = j['children'][0].quantitive_activity
            else:
                temp_dict['quantitive_activity'] = None
                temp_dict['quantitive_activity_initial'] = None
            temp_dict['qualitative_activity'] = j['children'][0].qualitative_activity
            temp_dict['quantitive_unit'] = j['children'][0].quantitive_unit
            temp_dict['quantitive_efficacy'] = j['children'][0].quantitive_efficacy
            temp_dict['efficacy_unit'] = j['children'][0].efficacy_unit
            temp_dict['quantitive_measure_type'] = j['children'][0].quantitive_measure_type
            temp_dict['efficacy_measure_type'] = j['children'][0].efficacy_measure_type
            temp_dict['t_coefficient'] = j['children'][0].bias_value
            temp_dict['t_coefficient_initial'] = j['children'][0].bias_value_initial
            temp_dict['bias_reference'] = j['children'][0].bias_reference
            temp_dict['emax_reference_ligand'] = j['children'][0].emax_ligand_reference
            temp_dict['ligand_function'] = j['children'][0].ligand_function
            temp_dict['ligand'] = j['main'].ligand
            if (temp_dict['quantitive_activity_initial'] and
                    temp_dict['quantitive_measure_type'] != "Effect at single point measurement"):
                temp_dict['quantitive_activity_initial'] = (-1) * math.log10(
                    temp_dict['quantitive_activity_initial'])
                temp_dict['quantitive_activity_initial'] = "{:.2F}".format(
                    Decimal(temp_dict['quantitive_activity_initial']))
            temp['ref_ligand_experiment'] = j['children'][0].emax_ligand_reference
            doubles.append(temp_dict)
            temp['assay'] = doubles
            send.append(temp)
        return send

    def combine_unique(self, data):
        '''
        combining tested assays and reference assays
        '''
        context = dict()
        print('data len', len(data))
        for j in data:
            name = str(j['publication'].id) + \
                '/'  + '/' + str(j['receptor'].id)
            temp_obj = list()
            if name in context:
                temp_obj = context[name]['assay']
            for i in j['assay']:
                temp_obj.append(i)
            context[name] = j
            context[name]['assay'] = temp_obj

        return context

    def process_referenced_assays(self, data):
        '''
        separate tested assays and reference assays
        '''
        assay_counter = 0
        for j in data.items():
            assays, reference = self.return_refenced_assays(j[1]['assay'])
            j[1].pop('assay')
            j[1]['assay_list'] = assays
            j[1]['reference_assays_list'] = reference
            assay_counter = assay_counter+len(j[1]['assay_list'])+len(j[1]['reference_assays_list'])
        print('assay_counter len', assay_counter)
        return data

    def return_refenced_assays(self, assays):
        # pylint: disable=no-member
        # no error
        main, reference = list(), list()
        for assay in assays:
            if assay['bias_reference'] != '':
                reference.append(assay)
            else:
                main.append(assay)
        sorted_main = sorted(main, key=lambda k: k['quantitive_activity']
                             if k['quantitive_activity'] else 999999, reverse=True)
        sorted_reference = reference
        if len(sorted_reference) == 0:
            self.get_reference_from_emax(assays)

        return sorted_main, sorted_reference

    def get_reference_from_emax(self, assays):
        reference_ligand = list()
        for i in assays:
            if i['emax_reference_ligand'] == i['ligand']:
                reference_ligand.append(i)
        return reference_ligand

    def separate_ligands(self, context):
        content = dict()
        for i in context.items():
            for assay in i[1]['assay_list']:
                name = str(i[1]['publication'].id) + \
                    '/'+ str(assay['ligand'].id) + '/' + str(i[1]['receptor'].id)
                if name in content:
                    content[name]['assay_list'].append(assay)
                else:
                    content[name] = dict()
                    content[name]['publication'] = i[1]['publication']
                    content[name]['ligand'] = assay['ligand']
                    content[name]['endogenous_ligand'] = i[1]['endogenous_ligand']
                    content[name]['receptor'] = i[1]['receptor']
                    content[name]['vendor_counter'] = i[1]['vendor_counter']
                    content[name]['authors'] = i[1]['authors']
                    content[name]['ref_ligand_experiment'] = i[1]['ref_ligand_experiment']
                    content[name]['article_quantity'] = i[1]['article_quantity']
                    content[name]['labs_quantity'] = i[1]['labs_quantity']
                    content[name]['assay_list'] = list()
                    content[name]['assay_list'].append(assay)
                    content[name]['reference_assays_list'] = i[1]['reference_assays_list']

        return content

# pylint: disable=C0301
    def limit_family_set(self, assay_list):
        # pylint: disable=no-member
        # no error
        families = list()
        proteins = set()
        for assay in assay_list:
            if assay['family'] not in proteins:
                proteins.add(assay['family'])
                families.append(assay)
            else:
                compare_val = next(item for item in families if item["family"] == assay['family'])
                try:
                    if assay['quantitive_activity'] < compare_val['quantitive_activity']:
                        families[:] = [d for d in families if d.get('family') != compare_val['family']]
                        families.append(assay)
                except:
                    continue
        return families

# pylint: disable=C0301
    def limit_family_set_subs(self, assay_list):
        families = list()
        proteins = set()
        for assay in assay_list:
            if assay['signalling_protein'] not in proteins:
                proteins.add(assay['signalling_protein'])
                families.append(assay)
            else:
                compare_val = next(item for item in families if item["signalling_protein"] == assay['signalling_protein'])
                # print('\n***dublicate', compare_val['signalling_protein'], compare_val['quantitive_activity'])
                try:
                    if assay['quantitive_activity'] < compare_val['quantitive_activity']:
                        families[:] = [d for d in families if d.get('signalling_protein') != compare_val['signalling_protein']]
                        families.append(assay)
                        # print('\n***swap', assay['signalling_protein'], assay['quantitive_activity'])
                except:
                    families.append(assay)
        print(proteins)
        return families

    def process_calculation(self, context):
        for i in context.items():
            test = dict()
            temp_obj = list()

            for j in i[1]['assay_list']:
                if j not in temp_obj:
                    temp_obj.append(j)
                else:
                    pass
            i[1]['assay_list'] = temp_obj


            test = sorted(i[1]['assay_list'], key=lambda k: k['quantitive_activity']
                          if k['quantitive_activity'] else 999999, reverse=False)

            for item in enumerate(test):
                item[1]['order_no'] = item[0]

            i[1]['biasdata'] = test
            i[1].pop('assay_list')
            # calculate log bias
            self.calc_bias_factor(i[1]['biasdata'], i[1]['reference_assays_list'])
            # recalculates lbf if it is negative
            i[1]['biasdata'] = self.validate_lbf(i)
            self.calc_potency_and_transduction(i[1]['biasdata'])

        return context

# pylint: disable=C0301
    def calc_bias_factor(self, biasdata, reference):
        most_reference = dict()
        most_potent = dict()
        for i in biasdata:
            if i['order_no'] == 0:
                most_potent = i
                most_reference = self.get_reference_assay(reference, most_potent)
                i['log_bias_factor'] = None

        for i in biasdata:
            if i['order_no'] != 0:
                temp_reference = self.get_reference_assay(reference, i)
                try:
                    if (i['quantitive_measure_type'].lower() == 'ec50' and temp_reference['quantitive_measure_type'].lower() == 'ec50' and
                            most_potent['quantitive_measure_type'].lower() == 'ec50' and most_reference['quantitive_measure_type'].lower() == 'ec50'):
                        a = 0
                        b = 0
                        c = 0
                        d = 0
                        a = math.log10(
                            most_potent['quantitive_efficacy'] / most_potent['quantitive_activity'])
                        b = math.log10(
                            most_reference['quantitive_efficacy'] / most_reference['quantitive_activity'])
                        c = math.log10(
                            i['quantitive_efficacy'] / i['quantitive_activity'])
                        d = math.log10(
                            temp_reference['quantitive_efficacy'] / temp_reference['quantitive_activity'])
                        temp_calculation = self.caclulate_bias_factor_variables(
                            a, b, c, d)
                        i['log_bias_factor'] = round(temp_calculation, 1)
                except:
                    if (i['quantitive_measure_type'].lower() == 'ic50' and temp_reference['quantitive_measure_type'].lower() == 'ic50'):
                        i['log_bias_factor'] = 'Only agonist in main pathway'
                    elif i['qualitative_activity'] == 'No activity':
                        i['log_bias_factor'] = "Full Bias"
                    elif i['qualitative_activity'] == 'Low activity':
                        i['log_bias_factor'] = "High Bias"

                    elif i['qualitative_activity'] == 'High activity':
                        i['log_bias_factor'] = "Low Bias"
                    else:
                        i['log_bias_factor'] = None

    def get_reference_assay(self, reference, assay):
        return_assay = dict()
        try:
            for i in reference:
                if i['signalling_protein'] == assay['signalling_protein']:
                    return_assay = i
        except:
            return return_assay
        return return_assay

    def caclulate_bias_factor_variables(self, a, b, c, d):
        '''
        calculations for log bias factor inputs
        '''
        lgb = 0
        try:
            lgb = (a - b) - (c - d)
        except:
            lgb = 0
        return lgb

    def calc_potency_and_transduction(self, biasdata):
        count = 0
        most_potent = dict()
        for i in biasdata:
            count += 1
            if i['order_no'] == 0:
                most_potent = i
        # T_factor -- bias factor
        for i in biasdata:
            if i['order_no'] > 0:
                if i['quantitive_measure_type'].lower() == 'ec50' or i['quantitive_measure_type'].lower() == 'ic50':

                    if i['quantitive_activity'] is not None and i['quantitive_activity'] != 0 and most_potent['quantitive_activity'] is not None:
                        i['potency'] = round(
                            i['quantitive_activity'] / most_potent['quantitive_activity'], 1)
                    elif i['quantitive_measure_type'].lower() == 'pec50' or i['quantitive_measure_type'].lower() == 'pic50':
                        i['potency'] = round(
                            most_potent['quantitive_activity'] - i['quantitive_activity'], 1)
                    else:
                        i['potency'] = None

                if i['t_coefficient'] is not None and most_potent['t_coefficient'] is not None:
                    i['t_factor'] = round(
                        most_potent['t_coefficient'] - i['t_coefficient'], 1)
                else:
                    i['t_factor'] = None

    def validate_lbf(self, i):
        for x in i[1]['biasdata']:
            if isinstance(x['log_bias_factor'], float):
                if x['log_bias_factor'] < 0.0:
                    j = next((item for item in i[1]['biasdata'] if item["order_no"] == 0), None)
                    x['order_no'], j['order_no'] = j['order_no'], x['order_no']
                    self.calc_bias_factor(i[1]['biasdata'], i[1]['reference_assays_list'])
                    self.validate_lbf(i)
                else:
                    return i[1]['biasdata']
        return i[1]['biasdata']

    def save_data_to_model(self, context, source):
        for i in context['data'].items():
            if self.fetch_experiment(i[1]['publication'], i[1]['ligand'], i[1]['receptor'], source) == False:
                primary, secondary = self.fetch_receptor_trunsducers(
                    i[1]['receptor'])

                experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],  # NOQA
                                                      ligand=i[1]['ligand'],
                                                      receptor=i[1]['receptor'],
                                                      source=source,
                                                      endogenous_ligand=i[1]['endogenous_ligand'],
                                                      vendor_quantity=i[1]['vendor_counter'],
                                                      reference_ligand=i[1]['ref_ligand_experiment'],
                                                      primary=primary,
                                                      secondary=secondary,
                                                      article_quantity=i[1]['article_quantity'],
                                                      labs_quantity=i[1]['labs'],
                                                      )
                experiment_entry.save()
                for ex in i[1]['biasdata']:
                    emax_ligand = ex['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     family=ex['family'],
                                                     order_no=ex['order_no'],
                                                     signalling_protein=ex['signalling_protein'],
                                                     cell_line=ex['cell_line'],
                                                     assay_type=ex['assay_type'],
                                                     assay_measure=ex['assay_measure_method'],
                                                     assay_time_resolved=ex['assay_time_resolved'],
                                                     ligand_function=ex['ligand_function'],
                                                     quantitive_measure_type=ex['quantitive_measure_type'],
                                                     quantitive_activity=ex['quantitive_activity'],
                                                     quantitive_activity_initial=ex['quantitive_activity_initial'],
                                                     quantitive_unit=ex['quantitive_unit'],
                                                     qualitative_activity=ex['qualitative_activity'],
                                                     quantitive_efficacy=ex['quantitive_efficacy'],
                                                     efficacy_measure_type=ex['efficacy_measure_type'],
                                                     efficacy_unit=ex['efficacy_unit'],
                                                     potency=ex['potency'],
                                                     t_coefficient=ex['t_coefficient'],
                                                     t_value=ex['t_coefficient_initial'],
                                                     t_factor=ex['t_factor'],
                                                     log_bias_factor=ex['log_bias_factor'],
                                                     emax_ligand_reference=emax_ligand
                                                     )
                    experiment_assay.save()
                for ex in i[1]['reference_assays_list']:
                    emax_ligand = ex['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     assay_description='reference_assay',
                                                     family=ex['family'],
                                                     order_no=ex['order_no'],
                                                     signalling_protein=ex['signalling_protein'],
                                                     cell_line=ex['cell_line'],
                                                     assay_type=ex['assay_type'],
                                                     assay_measure=ex['assay_measure_method'],
                                                     assay_time_resolved=ex['assay_time_resolved'],
                                                     ligand_function=ex['ligand_function'],
                                                     quantitive_measure_type=ex['quantitive_measure_type'],
                                                     quantitive_activity=ex['quantitive_activity'],
                                                     quantitive_activity_initial=ex['quantitive_activity_initial'],
                                                     quantitive_unit=ex['quantitive_unit'],
                                                     qualitative_activity=ex['qualitative_activity'],
                                                     quantitive_efficacy=ex['quantitive_efficacy'],
                                                     efficacy_measure_type=ex['efficacy_measure_type'],
                                                     efficacy_unit=ex['efficacy_unit'],
                                                     potency=ex['potency'],
                                                     t_coefficient=ex['t_coefficient'],
                                                     t_value=ex['t_coefficient_initial'],
                                                     t_factor=ex['t_factor'],
                                                     log_bias_factor=ex['log_bias_factor'],
                                                     emax_ligand_reference=emax_ligand
                                                     )
                    experiment_assay.save()

            else:
                pass

    def fetch_experiment(self, publication, ligand, receptor, source):
        '''
        fetch receptor with Protein model
        requires: protein id, source
        '''
        try:
            experiment = AnalyzedExperiment.objects.filter(
                publication=publication, ligand=ligand, receptor=receptor, source=source)
            experiment = experiment.get()
            return True
        except Exception:
            experiment = None
            return False

    def fetch_receptor_trunsducers(self, receptor):
        primary = set()
        temp = str()
        temp1 = str()
        secondary = set()
        try:
            gprotein = ProteinGProteinPair.objects.filter(protein=receptor)
            for x in gprotein:
                if x.transduction and x.transduction == 'primary':
                    primary.add(x.g_protein.name)
                elif x.transduction and x.transduction == 'secondary':
                    secondary.add(x.g_protein.name)
            for i in primary:
                temp += str(i) + str(', ')

            for i in secondary:
                temp1 += str(i) + str(', ')
            return temp, temp1
        except:
            print('receptor not found', receptor)
            return None, None

    def process_signalling_proteins(self, context):
        for i in context.items():
            i[1]['assay_list'] = self.limit_family_set(i[1]['assay_list'])
        return context

    def process_signalling_proteins_subs(self, context):
        for i in context.items():
            i[1]['assay_list'] = self.limit_family_set_subs(i[1]['assay_list'])
        return context

    def build_bias_data(self):
        print('Build bias data gproteins')
        context = dict()
        content = self.get_from_model()
        print('stage # 1 : Getting data finished, data points: ', len(content))
        content_with_children = self.process_data(content)
        print('stage # 2: Processing children in queryset finished', len(content_with_children))
        changed_data = self.queryset_to_dict(content_with_children)
        print('stage # 3: Converting queryset into dict finished', len(changed_data))
        send = self.combine_unique(changed_data)
        referenced_assay = self.process_referenced_assays(send)
        print('stage # 4: Separating reference assays is finished', len(referenced_assay))
        ligand_data = self.separate_ligands(referenced_assay)
        limit_family = self.process_signalling_proteins(ligand_data)
        print('stage # 5: Separate ligands', len(limit_family))
        calculated_assay = self.process_calculation(limit_family)
        print('stage # 6: Merging assays with same ligand/receptor/publication is finished')
        self.count_publications(calculated_assay)
        print('stage # 7: labs and publications counted')
        context.update({'data': calculated_assay})
        print('stage # 8: combining data into common dict is finished')
        # save dataset to model
        self.save_data_to_model(context, 'different_family')
        print('stage # 9: saving data to model is finished')

    def build_bias_data_subtypes(self):
        print('Build bias data gproteins')
        context = dict()
        content = self.get_from_model()
        print('stage # 1 : Getting data finished, data points: ', len(content))
        content_with_children = self.process_data(content)
        print('stage # 2: Processing children in queryset finished', len(content_with_children))
        changed_data = self.queryset_to_dict(content_with_children)
        print('stage # 3: Converting queryset into dict finished', len(changed_data))
        send = self.combine_unique(changed_data)
        referenced_assay = self.process_referenced_assays(send)
        print('stage # 4: Separating reference assays is finished', len(referenced_assay))
        ligand_data = self.separate_ligands(referenced_assay)
        limit_family = self.process_signalling_proteins_subs(ligand_data)
        print('stage # 5: Separate ligands', len(limit_family))
        calculated_assay = self.process_calculation(limit_family)
        print('stage # 6: Merging assays with same ligand/receptor/publication is finished')
        self.count_publications(calculated_assay)
        print('stage # 7: labs and publications counted')
        context.update({'data': calculated_assay})
        print('stage # 8: combining data into common dict is finished')
        # save dataset to model
        self.save_data_to_model(context, 'same_family')
        print('stage # 9: saving data to model is finished')

    def count_publications(self, context):
        temp = dict()
        for i in context.items():
            labs = list()
            i[1]['labs'] = 0
            labs.append(i[1]['publication'])
            lab_counter = 1
            for j in context.items():
                if j[1]['publication'] not in labs:
                    if set(i[1]['authors']) & set(j[1]['authors']):
                        lab_counter += 1
                        labs.append(j[1]['publication'])
                        i[1]['labs'] = lab_counter

            temp_obj = 1
            name = str(i[1]['ref_ligand_experiment']) + \
                '/' + str(i[1]['ligand']) + '/' + str(i[1]['receptor'])
            if name in temp:
                for assays in i[1]['biasdata']:
                    if assays['order_no'] > 0:
                        if assays['log_bias_factor'] != None and assays['log_bias_factor'] != '' or assays['t_factor'] != None and assays['t_factor'] != '':
                            temp_obj = temp[name] + 1

            temp[name] = temp_obj

        for i in context.items():
            temp_obj = 0
            name = str(i[1]['ref_ligand_experiment']) + \
                '/' + str(i[1]['ligand']) + '/' + str(i[1]['receptor'])
            if name in temp:
                i[1]['article_quantity'] = temp[name]
