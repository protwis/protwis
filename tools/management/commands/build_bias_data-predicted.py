from decimal import Decimal
import logging
import math
import pandas as pd
import os
from build.management.commands.base_build import Command as BaseBuild
from protein.models import ProteinGProteinPair
from ligand.models import BiasedExperiment, AnalyzedExperiment, AnalyzedAssay
from django.conf import settings


class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)
    structure_data_dir = os.sep.join(
        [settings.DATA_DIR, 'ligand_data', 'gproteins'])
    cell_structure_data_dir = os.sep.join(
        [settings.DATA_DIR, 'ligand_data', 'cell_line'])
    help = 'Reads bias data and imports it'

    gprot_cache = {}
    cell_cache = dict()

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
        self.logger.info('COMPLETED CREATING BIAS DATA')

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        # delete any existing structure data
        if options['purge']:
            try:
                print('Started purging bias data')
                Command.purge_bias_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
        print('CREATING BIAS DATA')
        print(options['filename'])
        self.build_bias_data()
        self.logger.info('COMPLETED CREATING BIAS DATA')

    @staticmethod
    def purge_bias_data():
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()


    @staticmethod
    def process_gproteins_excel():
        source_file_path = None
        filenames = os.listdir(Command.structure_data_dir)
        for source_file in filenames:
            source_file_path = os.sep.join(
                [Command.structure_data_dir, source_file]).replace('//', '/')
            print(source_file, source_file_path)
        df = pd.read_excel(source_file_path)
        Command.gprot_cache = df.set_index('UniProt').T.to_dict('dict')


    @staticmethod
    def process_cell_line_excel():
        source_file_path = None
        filenames = os.listdir(Command.cell_structure_data_dir)
        for source_file in filenames:
            source_file_path = os.sep.join(
                [Command.cell_structure_data_dir, source_file]).replace('//', '/')
            print(source_file, source_file_path)
        df = pd.read_excel(source_file_path)
        Command.cell_cache = df.set_index("Cell_line_name").T.to_dict('dict')

    def build_bias_data(self):
        print('prestage, process excell')
        Command.process_gproteins_excel()
        Command.process_cell_line_excel()
        print('Build bias data gproteins')
        context = dict()
        content = Command.get_data_from_model()
        # import pdb; pdb.set_trace()
        print('stage # 2 : Getting data finished, data points: ', len(content))
        content_with_children = Command.process_data(content)
        # import pdb; pdb.set_trace()
        print('stage # 3: Processing children in queryset finished',
              len(content_with_children))
        changed_data = self.queryset_to_dict(content_with_children)
        # import pdb; pdb.set_trace()
        print('stage # 4: Converting queryset into dict finished', len(changed_data))
        send = Command.combine_unique(changed_data)
        print('stage # 5: Selecting endogenous ligands finished')
        referenced_assay = self.process_referenced_assays(send)
        print('stage # 6: Separating reference assays is finished',
              len(referenced_assay))
        ligand_data = Command.separate_ligands(referenced_assay)
        print('stage # 7: Separate ligands finished')
        limit_family = self.process_signalling_proteins(ligand_data)
        print('stage # 8: process_signalling_proteins finished', len(limit_family))
        calculated_assay = Command.process_calculation(limit_family)
        print('stage # 9: Calucating finished')
        Command.count_publications(calculated_assay)
        print('stage # 10: labs and publications counted')
        context.update({'data': calculated_assay})
        print('stage # 11: combining data into common dict is finished')
        # save dataset to model
        Command.save_data_to_model(context, 'predicted_family')
        print('stage # 12: saving data to model is finished')

    @staticmethod
    def get_data_from_model():
        try:
            content = BiasedExperiment.objects.all().prefetch_related(
                'experiment_data', 'ligand', 'receptor', 'publication', 'publication__web_link', 'experiment_data__emax_ligand_reference',
            ).order_by('publication', 'receptor', 'ligand')
        except BiasedExperiment.DoesNotExist:
            content = None
        return content

    @staticmethod
    def process_data(content):
        rd = []
        counter = 0
        for instance in enumerate(content):
            temp_obj = []
            fin_obj = {}
            fin_obj['main'] = (instance[1])
            vendor_counter = 0
            vendors_quantity = None
            for i in instance[1].experiment_data_vendors.all():
                vendor_counter = vendor_counter + 1
                if not vendor_counter:
                    vendors_quantity = i

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

    @staticmethod
    def process_g_protein(protein, receptor):
        receptor_name = receptor.entry_name.split('_')[0].upper()
        if receptor_name in Command.gprot_cache:
            protein = Command.gprot_cache[receptor_name]["1'Gfam"]
        return protein

    @staticmethod
    def process_cell_line(cell_line):
        if cell_line in Command.cell_cache:
            _species = Command.cell_cache[cell_line]["Species"]
            _tissue = Command.cell_cache[cell_line]["Tissue/organ"]
        else:
            _species = cell_line
            _tissue = cell_line
        return _species, _tissue

    @staticmethod
    def queryset_to_dict(results):
        '''
        Merge bias experminet data with assay data
        '''
        send = list()
        for j in results:
            temp_dict = dict()
            temp = dict()
            temp['reference'] = list()
            temp['assay'] = dict()
            temp['ref_ligand_experiment'] = dict()
            doubles = []
            temp['publication'] = j['main'].publication
            temp['species'] = j['main'].receptor.species.common_name
            temp['ligand'] = j['main'].ligand
            temp['endogenous_ligand'] = j['main'].endogenous_ligand
            temp['auxiliary_protein'] = j['main'].auxiliary_protein
            temp['receptor'] = j['main'].receptor
            temp['receptor_isoform'] = j['main'].receptor_isoform
            temp['receptor_gtpo'] = j['main'].receptor_gtpo
            temp['vendor_counter'] = j['vendor_counter']
            temp['authors'] = j['authors']
            temp['article_quantity'] = 0
            temp['labs_quantity'] = 0
            if j['children']:
                temp_dict = Command.process_children_from_queryset(j, temp['receptor'])
                if temp_dict is not None:
                    doubles.append(temp_dict)
            temp['assay'] = doubles
            send.append(temp)
        return send

    @staticmethod
    def process_children_from_queryset(j, receptor):
        temp_dict = dict()
        temp_dict['assay_initial'] = j['children'][0]
        temp_dict['ligand_source_id'] = j['main'].ligand_source_id
        temp_dict['ligand_source_type'] = j['main'].ligand_source_type
        temp_dict['potency'] = None
        temp_dict['pathway_level'] = j['children'][0].pathway_level
        temp_dict['delta_relative_transduction_coef'] = None
        temp_dict['log_bias_factor'] = None
        temp_dict['delta_emax_ec50'] = None
        temp_dict['calculated_relative_tau'] = None
        temp_dict['order_no'] = 0
        temp_dict['endogenous_assay'] = dict() #shall be only one
        temp_dict['signalling_protein'] = j['children'][0].signalling_protein
        temp_dict['cell_line'] = j['children'][0].cell_line
        temp_dict['_tissue'], temp_dict['_species']  = Command.process_cell_line(temp_dict['cell_line'])
        temp_dict['family'] = j['children'][0].family
        if temp_dict['family'] == 'G protein' or temp_dict['family'] == 'Gq/11 or Gi/o':
            temp_dict['family'] = Command.process_g_protein(
                temp_dict['family'], receptor)
        if temp_dict['family'] == 'G protein' or temp_dict['family'] == 'Gq/11 or Gi/o':
            return None

        temp_dict['measured_biological_process'] = j['children'][0].measured_biological_process
        temp_dict['assay_type'] = j['children'][0].assay_type
        temp_dict['assay_time_resolved'] = j['children'][0].assay_time_resolved
        temp_dict['signal_detection_tecnique'] = j['children'][0].signal_detection_tecnique
        temp_dict['molecule_1'] = j['children'][0].molecule_1
        temp_dict['molecule_2'] = j['children'][0].molecule_2

        temp_dict['quantitive_activity'] = j['children'][0].quantitive_activity
        temp_dict['quantitive_activity_initial'] = j['children'][0].quantitive_activity

        temp_dict['qualitative_activity'] = j['children'][0].qualitative_activity
        temp_dict['quantitive_unit'] = j['children'][0].quantitive_unit
        temp_dict['quantitive_efficacy'] = j['children'][0].quantitive_efficacy
        temp_dict['efficacy_unit'] = j['children'][0].efficacy_unit
        temp_dict['quantitive_measure_type'] = j['children'][0].quantitive_measure_type
        temp_dict['efficacy_measure_type'] = j['children'][0].efficacy_measure_type
        temp_dict['transduction_coef'] = j['children'][0].transduction_coef
        temp_dict['relative_transduction_coef'] = j['children'][0].relative_transduction_coef

        temp_dict['bias_reference'] = j['children'][0].bias_reference
        temp_dict['emax_reference_ligand'] = j['children'][0].emax_ligand_reference
        temp_dict['ligand_function'] = j['children'][0].ligand_function
        temp_dict['ligand'] = j['main'].ligand
        temp_dict['quantitive_activity'], temp_dict['quantitive_activity_initial'] = Command.process_ec50_children_from_queryset(temp_dict)
        return temp_dict

    @staticmethod
    def process_ec50_children_from_queryset(temp_dict):
        try:
            temp_dict['quantitive_activity'] = float(temp_dict['quantitive_activity'])
        except:
            temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']
        if (temp_dict['quantitive_activity_initial'] and
                temp_dict['quantitive_measure_type'] != "Effect at single point measurement"):
            temp_dict['quantitive_activity_initial'] = (-1) * math.log10(
                temp_dict['quantitive_activity_initial'])
            temp_dict['quantitive_activity_initial'] = "{:.2F}".format(
                Decimal(temp_dict['quantitive_activity_initial']))
        return temp_dict['quantitive_activity'], temp_dict['quantitive_activity_initial']


    @staticmethod
    def combine_unique(data):
        '''
        combining tested assays and reference assays
        '''
        context = dict()
        for j in data:
            name = str(j['publication'].id) + \
                '/' + str(j['receptor'].id)
            temp_obj = list()
            if name in context:
                temp_obj = context[name]['assay']
            for i in j['assay']:
                temp_obj.append(i)
            context[name] = j
            context[name]['assay'] = temp_obj
        print("******len of experiments:", len(context), "******")
        return context

    @staticmethod
    def process_referenced_assays(data):
        '''
        separate tested assays and reference assays
        '''
        return_list = list()
        counter = 0
        for j in data.items():
            counter = counter+1
            return_dict = dict()
            assays, reference = Command.return_refenced_assays(j[1]['assay'])
            j[1]['assay_list'] = assays
            j[1]['reference_assays_list'] = reference
            j[1].pop('assay')
        return data

    _reference_assay_counter = 0
    @staticmethod
    def return_refenced_assays(assays):
        main, reference = list(), list()
        for assay in assays:
            # TODO: change to primary_Endogenous
            if (assay['bias_reference'] == 'Ref. and principal endo.' or
            assay['bias_reference'] == 'Endogenous' or
            assay['bias_reference'] == 'Principal endogenous' or
            assay['bias_reference'] ==  'Ref. and endo.'):
                if assay['quantitive_activity'] is not None:
                    reference.append(assay)
                    Command._reference_assay_counter = Command._reference_assay_counter+1
            else:
                main.append(assay)
        main = Command.fetch_endogenous_assay(main, reference)
        return main, reference

    @staticmethod
    def fetch_endogenous_assay(main, references):
        result_list = list()
        for assay in main:
            temp_reference_list = list()
            for reference in references:
                if assay['family'] == reference['family']:
                    # if assay['signalling_protein']:
                    if assay['signalling_protein'] == reference['signalling_protein']:
                        if assay['assay_type'] == reference['assay_type']:
                            if assay['cell_line'] == reference['cell_line']:
                                if assay['measured_biological_process'] == reference['measured_biological_process']:
                                    temp_reference_list.append(reference)
                    # else:
                    #     if assay['assay_type'] == reference['assay_type']:
                    #         if assay['cell_line'] == reference['cell_line']:
                    #             if assay['measured_biological_process'] == reference['measured_biological_process']:
                    #                 temp_reference_list.append(reference)

            if len(temp_reference_list)>0:
                if len(temp_reference_list)>1:
                    return_back_list = list()
                    final_end = None
                    for _reference_assay in temp_reference_list:
                        if _reference_assay['bias_reference'] == "Principal endogenous" or _reference_assay['bias_reference'] == "Ref. and principal endo.":
                            assay['endogenous_assay'] = _reference_assay
                            final_end = _reference_assay
                    if final_end is not None:
                        for _reference_assay in temp_reference_list:
                            if _reference_assay['bias_reference'] != "Principal endogenous" or _reference_assay['bias_reference'] != "Ref. and principal endo.":
                                _reference_assay['endogenous_assay'] = final_end
                                result_list.append(_reference_assay)
                else:
                    assay['endogenous_assay'] = temp_reference_list[0]

        return main

    @staticmethod
    def separate_ligands(context):
        content = dict()
        for i in context.items():
            if(len(i[1]['reference_assays_list'])) < 1:
                for assay in i[1]['assay_list']:
                    _pub_name = str(i[1]['publication'].id)
                    _ligand_name = str(assay['ligand'].id)
                    _receptor_name = str(i[1]['receptor'].id)
                    _receptor_iso_name = str(i[1]['receptor_isoform'])
                    _aux_prot_name = str(i[1]['auxiliary_protein'])
                    _tissue = assay['_tissue']
                    _species = assay['_species']
                    name = _pub_name+'/'+_ligand_name+'/'+_receptor_name+'/'+_receptor_iso_name+'/'+_aux_prot_name+'/'+_tissue+'/'+_species
                    if name in content:
                        content[name]['assay_list'].append(assay)

                    else:
                        content[name] = dict()
                        content[name]['assay_list'] = list()
                        content[name]['publication'] = i[1]['publication']
                        content[name]['ligand'] = assay['ligand']
                        content[name]['receptor_isoform']=i[1]['receptor_isoform']
                        content[name]['receptor_gtpo']=i[1]['receptor_gtpo']
                        content[name]['ligand_links'] = Command.get_external_ligand_ids(
                            content[name]['ligand'])
                        try:
                            content[name]['reference_ligand'] = i[1]['reference_assays_list'][0]['ligand']
                        except:
                            content[name]['reference_ligand'] = None
                        content[name]['auxiliary_protein'] = i[1]['auxiliary_protein']
                        # TODO: add external LigandStatistics
                        content[name]['endogenous_ligand'] = i[1]['endogenous_ligand']
                        content[name]['receptor'] = i[1]['receptor']
                        content[name]['vendor_counter'] = i[1]['vendor_counter']
                        content[name]['authors'] = i[1]['authors']
                        content[name]['article_quantity'] = i[1]['article_quantity']
                        content[name]['labs_quantity'] = i[1]['labs_quantity']
                        content[name]['assay_list'].append(assay)
                        content[name]['ligand_source_id'] = assay['ligand_source_id']
                        content[name]['ligand_source_type'] = assay['ligand_source_type']
        return content

    @staticmethod
    def get_external_ligand_ids(ligand):
        ligand_list = list()
        try:
            for i in ligand.properities.web_links.all():
                ligand_list.append(
                    {'name': i.web_resource.name, "link": i.index})
        except:
            ligand_list = list()
        return ligand_list

    @staticmethod
    def process_signalling_proteins(context):
        for i in context.items():

            i[1]['assay_list'] = Command.calculate_bias_factor_value(
                i[1]['assay_list'])

            i[1]['assay_list'] = Command.sort_assay_list(i[1]['assay_list'])

            i[1]['assay_list'] = Command.limit_family_set(i[1]['assay_list'])

            i[1]['assay_list'] = Command.order_assays(i[1]['assay_list'])
        return context

    @staticmethod
    def order_assays(assays):
        sorted_assay = sorted(assays, key=lambda k: k['delta_emax_ec50']
                                      if k['delta_emax_ec50'] else float(-1000), reverse=True)
        for item in enumerate(sorted_assay):
            item[1]['order_no'] = item[0]
        return assays

    @staticmethod
    def limit_family_set(assay_list):
        families = list()
        proteins = set()
        for assay in assay_list:
            if assay['family'] not in proteins:
                # if assay['family'] == 'Gq/11 or Gi/o':
                #     assay = self.get_rid_of_gprot(assay, families, proteins)
                if assay:
                    proteins.add(assay['family'])
                    families.append(assay)
            else:
                # if assay['family'] == 'Gq/11 or Gi/o':
                #     assay = self.get_rid_of_gprot(assay, families, proteins)
                if assay:
                    compare_val = next(
                        item for item in families if item["family"] == assay['family'])
                    try:
                        if assay['delta_emax_ec50'] > compare_val['delta_emax_ec50']:
                            families[:] = [d for d in families if d.get(
                                'family') != compare_val['family']]
                            families.append(assay)
                    except TypeError:
                        pass
        return families

    @staticmethod
    def sort_assay_list(i):
        return_assay = dict()
        return_assay = sorted(i, key=lambda k: k['delta_emax_ec50']
                              if k['delta_emax_ec50'] else float(-1000), reverse=True)
        return return_assay

    @staticmethod
    def calculate_bias_factor_value(sorted_assays):
        # TODO: pick
        for assay in sorted_assays:
            assay['delta_emax_ec50'] = Command.calc_delta_emax_ec50(assay)
        return sorted_assays

    @staticmethod
    def calc_delta_emax_ec50(assay):
        result = None
        try:
            assay_a = assay['quantitive_activity']
            assay_b = assay['quantitive_efficacy']
            result = math.log10(assay_b / assay_a)
        except:
            result = None
        return result

    @staticmethod
    def process_calculation(context):
        list_to_remove = list()
        for i in context.items():
            if len(i[1]['assay_list'])>1:
                i[1]['biasdata'] = i[1]['assay_list']
                i[1].pop('assay_list')
                # calculate log bias
                Command.calc_bias_factor(i[1]['biasdata'])
            else:
                list_to_remove.append(i[0])
        for experiment in list_to_remove:
            context.pop(experiment)
        return context

    @staticmethod
    def calc_bias_factor(biasdata):
        most_potent = dict()
        for i in biasdata:
            if i['order_no'] == 0:
                most_potent = i

        for i in biasdata:
            if i['order_no'] != 0:
                i['log_bias_factor'] = Command.lbf_process_qualitative_data(i)
                if i['log_bias_factor'] == None:
                    i['log_bias_factor'] = Command.lbf_process_ic50(i)
                if i['log_bias_factor'] == None:
                    try:
                        i['log_bias_factor'] = round(most_potent['delta_emax_ec50'] - i['delta_emax_ec50'], 1)
                    except:
                        i['log_bias_factor'] = None

    @staticmethod
    def lbf_process_qualitative_data(i):
        return_message = None
        try:
            if i['qualitative_activity'] == 'No activity':
                return_message = "Full Bias"
            elif i['qualitative_activity'] == 'Low activity':
                return_message = "High Bias"
            elif i['qualitative_activity'] == 'High activity':
                return_message = "Low Bias"
            elif i['qualitative_activity'] == 'Inverse agonism/antagonism':
                return_message = "Full Bias"
        except:
            return_message = None
        return return_message


    @staticmethod
    def lbf_process_ic50(i):
        return_message = None
        try:
            if (i['quantitive_measure_type'].lower() == 'ic50'):
                return_message = 'Only agonist in main pathway'
        except:
            return_message = None
        return return_message

    @staticmethod
    def count_publications(context):
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
            name = str(i[1]['endogenous_ligand']) + \
                '/' + str(i[1]['ligand'])+'/'+str(i[1]['receptor'])
            if name in temp:
                for assays in i[1]['biasdata']:
                    if assays['order_no'] > 0:
                        if assays['log_bias_factor'] != None and assays['log_bias_factor'] != '' or assays['delta_relative_transduction_coef'] != None and assays['delta_relative_transduction_coef'] != '':
                            temp_obj = temp[name] + 1

            temp[name] = temp_obj

        for i in context.items():
            temp_obj = 0
            name = str(i[1]['endogenous_ligand']) + \
                '/' + str(i[1]['ligand'])+'/'+str(i[1]['receptor'])
            if name in temp:
                i[1]['article_quantity'] = temp[name]

    @staticmethod
    def save_data_to_model(context, source):
        for i in context['data'].items():
            if len(i[1]['biasdata']) > 1:
                experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],
                                                      ligand=i[1]['ligand'],
                                                      external_ligand_ids=i[1]['ligand_links'],
                                                      receptor=i[1]['receptor'],
                                                      source=source,
                                                      receptor_isoform=i[1]['receptor_isoform'],
                                                      receptor_gtpo=i[1]['receptor_gtpo'],
                                                      endogenous_ligand=i[1]['endogenous_ligand'],
                                                      vendor_quantity=i[1]['vendor_counter'],
                                                      reference_ligand=i[1]['reference_ligand'],
                                                      article_quantity=i[1]['article_quantity'],
                                                      labs_quantity=i[1]['labs'],
                                                      ligand_source_id=i[1]['ligand_source_id'],
                                                      ligand_source_type=i[1]['ligand_source_type']
                                                      )
                experiment_entry.save()
                for ex in i[1]['biasdata']:
                    emax_ligand = ex['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     assay_description='predicted_tested_assays',
                                                     family=ex['family'],
                                                     order_no=ex['order_no'],
                                                     signalling_protein=ex['signalling_protein'],
                                                     cell_line=ex['cell_line'],
                                                     assay_type=ex['assay_type'],
                                                     pathway_level=ex['pathway_level'],
                                                     reference_assay_initial=None,
                                                     molecule_1=ex['molecule_1'],
                                                     molecule_2=ex['molecule_2'],
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
                                                     relative_transduction_coef=ex['relative_transduction_coef'],
                                                     transduction_coef=ex['transduction_coef'],
                                                     delta_relative_transduction_coef=ex['delta_relative_transduction_coef'],
                                                     log_bias_factor=ex['log_bias_factor'],
                                                     delta_emax_ec50=ex['delta_emax_ec50'],
                                                     effector_family=ex['family'],
                                                     measured_biological_process=ex['measured_biological_process'],
                                                     signal_detection_tecnique=ex['signal_detection_tecnique'],
                                                     emax_ligand_reference=emax_ligand
                                                     )
                    experiment_assay.save()
