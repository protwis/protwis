from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_or_create_ligand
from protein.models import Protein
from ligand.models import BiasedPathwaysAssay, BiasedPathways
from common.models import Publication
from common.tools import test_model_updates
from django.conf import settings
import logging
import os
import xlrd
import django.apps

class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    help = 'Reads bias data and imports it'
    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'bias_data', 'pathways'])
    publication_cache = {}
    ligand_cache = {}
    data_all = []
    mol_types = {"PubChem CID": "pubchem",
                "ChEMBL Compound ID": "chembl_ligand",
                "IUPHAR/BPS Guide to pharmacology": "gtoplig",
                "SMILES" : "smiles",
                "FASTA sequence (peptide)": "sequence",
                "UniProt Entry Code (peptide)" : "uniprot"}

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
                self.purge_bias_data()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        print('CREATING BIAS PATHWAYS DATA')
        self.prepare_all_data(options['filename'])
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('COMPLETED CREATING BIAS')

    def purge_bias_data(self):
        BiasedPathwaysAssay.objects.all().delete()
        BiasedPathways.objects.all().delete()

    def loaddatafromexcel(self, excelpath):
        """Reads excel file (require specific excel sheet)."""
        num_rows = 0
        try:
            workbook = xlrd.open_workbook(excelpath)
            worksheets = workbook.sheet_names()
            worksheets.remove("Dropdown") # remove annotation info sheet
            temp = []

            for worksheet_name in worksheets:
                worksheet = workbook.sheet_by_name(worksheet_name)
                num_rows = worksheet.nrows - 1
                num_cells = worksheet.ncols - 1
                curr_row = 0  # skip first, otherwise -1
                while curr_row < num_rows:
                    curr_row += 1
                    #row = worksheet.row(curr_row)
                    curr_cell = -1
                    temprow = []
                    while curr_cell < num_cells:
                        curr_cell += 1
                        cell_value = worksheet.cell_value(curr_row, curr_cell)
                        #cell_type = worksheet.cell_type(curr_row, curr_cell)

                        # fix wrong spaced cells
                        if cell_value == " ":
                            cell_value = ""
                        temprow.append(cell_value)
                    temp.append(temprow)
                    # if curr_row>10: break
            return temp
        except:
            self.logger.info(
                "The error appeared during reading the excel" + num_rows)

    def analyse_rows(self, rows):
        """
        Reads excel rows one by one.
        Fetch data to models.
        Saves to DB.
        """
        # Analyse the rows from excel and assign the right headers
        temp = []
        for i, r in enumerate(rows, 1):
            d = {}
            if r[4] != '':  # checks if the receptor field exists
                # try:
                d['submitting_group'] = r[0]
                # doi
                d['reference'] = r[1]
                # protein
                d['receptor'] = r[4].lower()
                d['signalling_protein'] = r[5].lower().strip()
                # ligand
                d['ligand_name'] = r[6]
                d['ligand_type'] = r[7]
                d['ligand_id'] = r[8]
                #experiment
                d['experiment_distinction'] = r[9]
                d['experiment_system'] = r[10]
                d['experiment_method'] = r[11]
                # pathway
                d['pathway_detail'] = r[12]
                d['pathway_summary'] = r[13]
                d['pathway_outcome'] = r[14]
                #Therapeutic
                d['effect_type'] = r[15]
                d['relevance'] = r[16]

                if not isinstance(d['ligand_id'], str):
                    d['ligand_id'] = int(d['ligand_id'])

                # fetch publicaition
                pub = self.fetch_publication(d['reference'])

                # fetch main ligand
                ligand = None
                chembl = None
                if d['ligand_name'] is not None and d['ligand_name'] != "":
                    ids = {}
                    if d['ligand_type'] in self.mol_types:
                        ids = {self.mol_types[d['ligand_type']]: d['ligand_id']}
                    ligand = get_or_create_ligand(d['ligand_name'], ids)
                    chembl = self.fetch_chembl(ligand)

                # fetch protein
                protein = self.fetch_protein(d['receptor'])
                if protein == None:
                    continue

                ## TODO:  check if it was already uploaded
                experiment_entry = BiasedPathways(submission_author=d['submitting_group'],
                                                    publication=pub,
                                                    ligand=ligand,
                                                    receptor=protein,
                                                    chembl = chembl,
                                                    relevance = d['relevance'],
                                                    effect_type = d['effect_type'],
                                                    signalling_protein = d['signalling_protein']
                                                    )
                experiment_entry.save()

                experiment_assay = BiasedPathwaysAssay(biased_pathway=experiment_entry,
                                                  pathway_outcome_high = d['pathway_outcome'],
                                                  pathway_outcome_summary = d['pathway_summary'],
                                                  pathway_outcome_detail  = d['pathway_detail'],
                                                  experiment_pathway_distinction = d['experiment_distinction'],
                                                  experiment_system = d['experiment_system'],
                                                  experiment_outcome_method= d['experiment_method']
                                                   )

                experiment_assay.save()

                # except Exception as msg:
                #     print(d['source_file'], msg)
                #     continue

            temp.append(d)
        return temp

    def fetch_chembl(self,ligand):
        links = ligand.ids.filter(web_resource__slug="chembl_ligand")
        if links.exists() > 0:
            return links.first().index
        else:
            return None

    def fetch_protein(self, protein_from_excel):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        test = None
        if Protein.objects.filter(entry_name=protein_from_excel).exists():
            return Protein.objects.get(entry_name=protein_from_excel)
        elif Protein.objects.filter(web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot').exists():
            return Protein.objects.get(web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot')
        else:
            return test

    def fetch_publication(self, publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        try:
            float(publication_doi)
            publication_doi = str(int(publication_doi))
        except ValueError:
            pass

        if publication_doi.isdigit():  # assume pubmed
            pub_type = 'pubmed'
        else:  # assume doi
            pub_type = 'doi'

        if publication_doi not in self.publication_cache:
            pub = False
            if pub_type == 'doi':
                pub = Publication.get_or_create_from_doi(publication_doi)
            elif pub_type == 'pubmed':
                pub = Publication.get_or_create_from_pubmed(publication_doi)

            if not pub:
                self.mylog.debug(
                    "publication fetching error | module: fetch_publication. Row # is : " + str(publication_doi) + ' ' + pub_type)

            self.publication_cache[publication_doi] = pub
        else:
            pub = self.publication_cache[publication_doi]

        return pub

    def prepare_all_data(self, filenames):
        if not filenames:
            filenames = os.listdir(self.structure_data_dir)
        for source_file in filenames:
            # print("source_file " + str(source_file))
            source_file_path = os.sep.join([self.structure_data_dir, source_file]).replace('//', '/')
            # print("source_file_path " + str(source_file_path))
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                print('Reading file {}'.format(source_file_path))
                # read the yaml file
                rows = []
                if source_file[-4:] == 'xlsx' or source_file[-3:] == 'xls':
                    if "~$" in source_file:
                        # ignore open excel files
                        continue
                    rows = self.loaddatafromexcel(source_file_path)
                    rows = self.analyse_rows(rows)
                else:
                    self.mylog.debug('unknown format'.source_file)
                    continue

                self.data_all += rows
        print(len(self.data_all), " total data points")
        print("Finished")
