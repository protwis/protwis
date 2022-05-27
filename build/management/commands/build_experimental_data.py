from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *
from django.conf import settings
from django.db.models import Prefetch
from django.utils.text import slugify

from common.models import WebLink, WebResource, Publication
from ligand.models import Ligand, LigandType, ExperimentalData, LigandVendors, LigandVendorLink #AssayExperiment
from protein.models import Protein

import math
import statistics
import datamol as dm
import datetime
import pandas as pd
import requests
import xmltodict
import urllib.parse
import urllib.request

class Command(BaseBuild):
    help = "Build ChEMBL, PDSP and GtoP bioactivities data"
    bulk_size = 50000

    def add_arguments(self, parser):
        parser.add_argument("--test_run",
                            action="store_true",
                            help="Skip this during a test run",
                            default=False)
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing ligand records')

    def handle(self, *args, **options):
        if options["test_run"]:
            print("Skipping in test run")
            return

        if options['purge']:
            # delete any existing ChEMBL bioactivity data
            # For deleting the ligands - purse all ligand data using the GtP ligand build
            print("Started purging bioactivity data")
            self.purge_experimental_data()
            print("Ended purging bioactivity data")

        # Parse ChEMBL bioactivity data
        print("\n\nStarted building ChEMBL bioactivities")
        self.build_chembl_bioactivities()
        print("Ended building ChEMBL bioactivities")

        #Fetching all the Guide to Pharmacology data
        print("\n\nStarted parsing Guide to Pharmacology bioactivities data")
        gtp_uniprot_link = "https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv"
        gtp_uniprot = pd.read_csv(gtp_uniprot_link, dtype=str, header=1)
        gtp_complete_ligands_link = "https://www.guidetopharmacology.org/DATA/ligands.csv"
        gtp_complete_ligands = pd.read_csv(gtp_complete_ligands_link, dtype=str, header=1)
        gtp_ligand_mapping_link = "https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv"
        gtp_ligand_mapping = pd.read_csv(gtp_ligand_mapping_link, dtype=str, header=1)
        gtp_interactions_link = "https://www.guidetopharmacology.org/DATA/interactions.csv"
        gtp_interactions = pd.read_csv(gtp_interactions_link, dtype=str, header=1)
        gtp_peptides_link = "https://www.guidetopharmacology.org/DATA/peptides.csv"
        gtp_peptides = pd.read_csv(gtp_peptides_link, dtype=str, header=1)
        # This gets all the info of the ligand and the interaction with the target
        iuphar_ids = self.compare_proteins(gtp_uniprot)
        bioactivity_ligands_ids = self.obtain_ligands(gtp_interactions, iuphar_ids, ['target_id','ligand_id'])
        #Now I have all the data I need
        bioactivity_data_gtp = self.get_bioligands_data(bioactivity_ligands_ids, gtp_complete_ligands, gtp_ligand_mapping, gtp_interactions)
        #Assess the assay type given info from affinity units and assay comments
        bioactivity_data_gtp = self.classify_assay(bioactivity_data_gtp)
        bioactivity_data_gtp.fillna('None', inplace=True)
        print("Ended parsing Guide to Pharmacology bioactivities data")

        # Building GtP bioactivity data
        print("\n\nStarted building Guide to Pharmacology bioactivities")
        build_gtp_bioactivities(bioactivity_data_gtp)
        print("Ended building Guide to Pharmacology bioactivities")

        # Building PDSP KiDatabase bioactivity data
        print("\n\nStarted building PDSP KiDatabase bioactivities")
        build_kidatabase_bioactivities() #14,562
        print("Ended building PDSP KiDatabase bioactivities")

    @staticmethod
    def purge_experimental_data():
        delete_experimental = ExperimentalData.objects.all() #New Model Biased Data
        delete_experimental.delete()

    @staticmethod
    def build_chembl_bioactivities():
        print("\n===============\n#1 Reading ChEMBL bioacitivity data")
        bioactivity_input_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_bioactivity_data.csv.gz"])
        bioactivity_data = pd.read_csv(bioactivity_input_file, dtype=str)
        url_template = 'https://www.ebi.ac.uk/chembl/api/data/document?document_chembl_id={}'
        bioactivity_data.fillna('None', inplace=True)
        bio_entries = len(bioactivity_data)
        print("Found", bio_entries, "bioactivities", datetime.datetime.now())

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        #ids = list(bioactivity_data["parent_molecule_chembl_id"].unique())  # not filtering is way faster
        ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Building ChEMBL proteins cache", datetime.datetime.now())
        # NOTE => might need to switch to Accession as the Entry name changes more frequently
        # If so, keep isoform notations in mind
        names = list(bioactivity_data["Entry name"].unique())
        proteins = list(Protein.objects.filter(entry_name__in=names).values_list("pk", "entry_name"))
        prot_dict = {prot_entry[1]: prot_entry[0] for prot_entry in proteins}

        print("\n#4 Building ChEMBL bioactivity entries", datetime.datetime.now())
        bioacts = []
        pub_links = []
        for index, row in bioactivity_data.iterrows():
            try:
                if (row["parent_molecule_chembl_id"] in lig_dict.keys()) and (row["Entry name"] in prot_dict.keys()):
                    bioacts.append(ExperimentalData())
                    bioacts[-1].ligand_id = lig_dict[row["parent_molecule_chembl_id"]]
                    bioacts[-1].protein_id = prot_dict[row["Entry name"]]
                    bioacts[-1].assay_type = row["assay_type"]
                    bioacts[-1].assay_description = row["assay_description"]
                    bioacts[-1].standard_activity_value = row["standard_value"]
                    bioacts[-1].p_activity_value = row["pchembl_value"]
                    bioacts[-1].p_activity_ranges = None
                    bioacts[-1].standard_relation = row["standard_relation"]
                    bioacts[-1].value_type = row["standard_type"]
                    bioacts[-1].document_chembl_id = row["document_chembl_id"]
                    bioacts[-1].source = 'ChEMBL'

                    response = requests.get(url_template.format(row["document_chembl_id"))
                    try:
                        data = xmltodict.parse(response.content)
                        doi = data['response']['documents']['document']['doi']
                        if doi is not None:
                            publication = fetch_publication(doi)
                            if ((len(bioacts)-1), publication.id) not in pub_links:
                                pub_links.append((len(bioacts)-1), publication.id)
                    except:
                        pass

                    # BULK insert every X entries or last entry
                    if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                        ExperimentalData.objects.bulk_create(bioacts)
                        data_pub = ExperimentalData.publication.through
                        data_pub.objects.bulk_create([
                                            data_pub(experimentaldata_id=data_id, publication_id=pub_id)
                                            for (data_id, pub_id) in pub_links])
                        print("Inserted", index, "out of", bio_entries, "bioactivities")
                        bioacts = []
                        pub_links = []
            except KeyError:
                continue

    @staticmethod
    def uniprot_mapper(protein, organism):
        organism_dict = {'PIG': 'sus_scrofa','RAT': 'rattus_norvegicus','HUMAN': 'homo_sapiens','MOUSE': 'mus_musculus',
                         'CANINE': 'canis_lupus_familiaris','BOVINE': 'bos_taurus','CALF': 'bos_taurus','COW': 'bos_taurus',
                         'GUINEA PIG': 'cavia_porcellus', 'CAT': 'felis_catus','NEONATAL RAT': 'rattus_norvegicus',
                         '? HUMAN': 'homo_sapiens','OPOSSUM': 'didelphis_marsupialis', 'Rat 6B': 'rattus_norvegicus',
                         'HUMAN M3': 'homo_sapiens','HUMAN M4': 'homo_sapiens', 'Chick': 'gallus_gallus',
                         'Frog': 'pseudis_balbodactyla','Newborn rats': 'rattus_norvegicus', 'Beef': 'bos_taurus',
                         'Sheep': 'ovis_aries','OX': 'bos_taurus','Dog': 'canis_lupus_familiaris',
                         'Rhesus': 'macaca_mulatta','Monkey': 'macaca_mulatta','PIGLET': 'sus_scrofa',
                         'Rat Y861': 'rattus_norvegicus','Zebra Finch': 'taeniopygia_guttata','Chicken': 'gallus_gallus',
                         'MICE': 'mus_musculus','Rhesus Monkey': 'macaca_mulatta', 'Zebrafish': 'danio_rerio'}
        if organism in organism_dict.keys():
            query = 'gene_exact:{0} AND organism:{1}'.format(protein.lower(), organism_dict[organism])
        else:
            query = 'gene_exact:{}'.format(protein.lower())

        url = 'https://www.uniprot.org/uniprot/'
        params = {
            'query': query,
            'format': 'tab',
            'columns': 'entry_name'
        }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        try:
            converted = urllib.request.urlopen(req).read().decode('utf-8').split('\n')[1].lower()
        except IndexError:
            converted = None

        return converted


    @staticmethod
    def get_bioligands_data(ligands, complete_ligands, ligand_mapping, ligand_interactions):
        full_info = ['Ligand id', 'Name','Species','Type','Approved','Withdrawn','Labelled','Radioactive', 'PubChem SID', 'PubChem CID',
                     'UniProt id','IUPAC name', 'INN', 'Synonyms','SMILES','InChIKey','InChI','GtoImmuPdb','GtoMPdb']
        bioactivity_info = ['ligand_id', 'action','target', 'target_id', 'target_species',
                            'target_gene_symbol', 'target_uniprot', 'original_affinity_relation',
                            'action_comment', 'selectivity', 'primary_target',
                            'concentration_range', 'affinity_units', 'affinity_high',
                            'affinity_median', 'affinity_low', 'assay_description', 'pubmed_id']
        weblinks = ['Ligand id', 'ChEMBl ID','Chebi ID','CAS','DrugBank ID','Drug Central ID']

        ligand_data = complete_ligands.loc[complete_ligands['Ligand id'].isin(ligands), full_info]
        ligand_weblinks = ligand_mapping.loc[ligand_mapping['Ligand id'].isin(ligands), weblinks]
        bioactivity_data = ligand_interactions.loc[ligand_interactions['ligand_id'].isin(ligands), bioactivity_info]
        bioactivity_data = bioactivity_data.rename(columns={'ligand_id': 'Ligand id'})

        ligand_complete = ligand_data.merge(ligand_weblinks, on="Ligand id")
        ligand_complete = ligand_complete.merge(bioactivity_data, on="Ligand id")
        ligand_complete = ligand_complete.rename(columns={"Ligand id": "Ligand ID"})

        return ligand_complete

    @staticmethod
    def obtain_ligands(data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(interactions_targets).intersection(set(compare_set))
        ligands = list(set(data.loc[data[labels[0]].isin(targets_with_ligands), labels[1]]))
        ligands = [x for x in ligands if x == x] #remove nan
        return ligands

    @staticmethod
    def compare_proteins(gtp_data):
        #Probably the files have been changed:
        #Now "uniprot_id" is "UniProtKB ID"
        #and "iuphar_id" is "GtoPdb IUPHAR ID"
        gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name','accession')
        # entries = gtp_data.loc[gtp_data['uniprot_id'].isin([protein[1].split("-")[0] for protein in gpcrdb_proteins]), ['uniprot_id', 'iuphar_id']]
        entries = gtp_data.loc[gtp_data['UniProtKB ID'].isin([protein[1].split("-")[0] for protein in gpcrdb_proteins]), ['UniProtKB ID', 'GtoPdb IUPHAR ID']]
        # return list(entries['iuphar_id'].unique())
        return list(entries['GtoPdb IUPHAR ID'].unique())

    @staticmethod
    def classify_assay(biodata):
        #starting to assess assay type
        biodata['assay_type'] = 'U'
        biodata['assay_description'].fillna('Unclassified',inplace=True)
        # find Binding assays (displacement, affinity, ) KI and KD are always binding
        biodata.loc[biodata['affinity_units'].isin(['pKi', 'pKd']), 'assay_type'] = 'B'
        # EC50, pKB and pA2 are always Functional
        biodata.loc[biodata['affinity_units'].isin(['pKB', 'pEC50', 'pA2']), 'assay_type'] = 'F'
        # IC50 can be both B: displacement, affinity, binding, radioligand else is F
        binding_words = ['isplacement', 'ffinity', 'inding', 'adioligand']
        biodata.loc[(biodata['affinity_units'] == 'pIC50') & (biodata['assay_description'].str.contains('|'.join(binding_words))), 'assay_type'] = 'B'
        biodata.loc[(biodata['affinity_units'] == 'pIC50') & ~(biodata['assay_description'].str.contains('|'.join(binding_words))), 'assay_type'] = 'F'
        # find Unclassified assayas
        biodata.loc[biodata['assay_description'].str.contains('Unclassified'), 'assay_type'] = 'U'

        return biodata

    @staticmethod
    def build_gtp_bioactivities(gtp_biodata):
        print("# Start parsing the GTP Dataframe")
        for index, row in gtp_biodata.iterrows():
            receptor = Command.fetch_protein(row['target_id'], row['target_species'])
            # TODO Handle multiple matches (uniprot filter?)
            ligand = get_ligand_by_id("gtoplig", row['Ligand ID'])

            try:
                low_value = "{:.2f}".format(float(row['affinity_low']))
            except ValueError:
                low_value = 'None'

            try:
                high_value = "{:.2f}".format(float(row['affinity_high']))
            except ValueError:
                high_value =  'None'

            if row['affinity_median'] != 'None':
                activity_value = "{:.2f}".format(float(row['affinity_median']))
            elif row['affinity_high'] != 'None':
                if row['affinity_low'] != 'None':
                    activity_value = "{:.2f}".format(statistics.mean([float(row['affinity_high']),float(row['affinity_low'])]))
                else:
                    activity_value = "{:.2f}".format(float(row['affinity_high']))
            else:
                activity_value = 'None'

            ranges = '|'.join([low_value, activity_value, high_value])

            #Adding publications from the PMIDs section
            try:
                pmids = row['pubmed_id'].split('|')
            except AttributeError:
                pmids = None

            if (receptor is not None) and (ligand is not None):
            #last step because it requires multiple uploads in case we have multiple species
                    gtp_data = ExperimentalData(
                                ligand = ligand,
                                protein = receptor,
                                assay_type = row['assay_type'],
                                assay_description = row['assay_description'],
                                standard_activity_value = None,
                                p_activity_value = activity_value,
                                p_activity_ranges = ranges,
                                standard_relation = row['original_affinity_relation'],
                                value_type = row['affinity_units'],
                                source = 'Guide to Pharmacology',
                                document_chembl_id = None,
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication= None
            else:
                print("SKIPPING", row["Ligand ID"], row['target_id'])

    @staticmethod
    def fetch_publication(publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid

        """
        if ("ISBN" in publication_doi) or (int(publication_doi) == 0):
            return None

        try:
            float(publication_doi)
            publication_doi = str(int(publication_doi))
        except ValueError:
            pass

        if publication_doi.isdigit():  # assume pubmed
            pub_type = 'pubmed'
        else:  # assume doi
            pub_type = 'doi'

        try:
            wl = WebLink.objects.get(index=publication_doi, web_resource__slug=pub_type)
        except WebLink.DoesNotExist:
            try:
                wl = WebLink.objects.create(index=publication_doi, web_resource=WebResource.objects.get(slug=pub_type))
            except IntegrityError:
                wl = WebLink.objects.get(index=publication_doi, web_resource__slug=pub_type)

        try:
            pub = Publication.objects.get(web_link=wl)
        except Publication.DoesNotExist:
            pub = Publication()
            try:
                pub.web_link = wl
                pub.save()
            except IntegrityError:
                pub = Publication.objects.get(web_link=wl)

            if pub_type == 'doi':
                pub.update_from_doi(doi=publication_doi)
            elif pub_type == 'pubmed':
                pub.update_from_pubmed_data(index=publication_doi)
            try:
                pub.save()
            except:
                self.mylog.debug(
                    "publication fetching error | module: fetch_publication. Row # is : " + str(publication_doi) + ' ' + pub_type)
                # if something off with publication, skip.
        return pub

    @staticmethod
    def fetch_protein(target, species):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            if species == None or species == "None":
                # Sorting by species => human first, otherwise next species in line
                # TODO => potentially capture all species with GtP ID
                return Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop').order_by("species_id").first()
            else:
                prots = Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop', species__common_name__iexact=species)
                if prots.count() > 0:
                    return prots.first()
                else:
                    receptor_fam = list(Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop').values_list("family_id", flat = True))[0]
                    return Protein.objects.get(family_id=receptor_fam, species__common_name__iexact=species)
        except:
            return None

    @staticmethod
    def fetch_protein_from_kidata(protein_input):
        """
        fetch receptor with Protein model
        requires: protein entry name
        """
        try:
            test = None
            if Protein.objects.filter(entry_name=protein_input):
                protein = Protein.objects.filter(entry_name=protein_input)
                test = protein.get()
            elif Protein.objects.filter(web_links__index=protein_input, web_links__web_resource__slug='uniprot'):
                protein1 = Protein.objects.filter(web_links__index=protein_input, web_links__web_resource__slug='uniprot')
                test = protein1[0]
            return test
        except:
            return None

    @staticmethod
    def build_kidatabase_bioactivities():
        protein_names = {}
        ligand_cache = {}
        print("\n===============\n#1 Reading PDSP bioacitivity data")
        bioactivity_kidatabase_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "KiDatabase.csv"])
        bioactivity_kidata = pd.read_csv(bioactivity_kidatabase_file, dtype=str)
        #Keeping data that has either SMILES info OR CAS info
        #CAS number can be translated into pubchem CID
        bioactivity_data_filtered = bioactivity_kidata.loc[(~bioactivity_kidata['SMILES'].isnull()) | (~bioactivity_kidata['CAS'].isnull())]
        bioactivity_data_filtered = bioactivity_data_filtered.loc[(~bioactivity_data_filtered['Unigene'].isnull())]
        bioactivity_data_filtered.fillna('None', inplace=True)
        bio_entries = len(bioactivity_data_filtered)
        print("\n===============\n#2 Start parsing PDSP data")
        bioacts = []
        for index, row in bioactivity_data_filtered.iterrows():
            ids = {}
            label = '_'.join([row['Unigene'], row['species']])
            if label not in protein_names.keys():
                protein = Command.uniprot_mapper(row['Unigene'], row['species'])
                if protein is not None:
                    protein_names[label] = protein
                else:
                    continue
            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
            if row['CAS'] != 'None':
                ids['CAS'] = row['CAS']
            if row[' Ligand Name'] not in ligand_cache.keys():
                ligand = get_or_create_ligandrow[' Ligand Name'], ids)
                ligand_cache[row[' Ligand Name']] = ligand
            receptor = Command.fetch_protein_from_kidata(protein_names[label])
            if (receptor is not None) and (ligand_cache[row[' Ligand Name']] is not None):
                bioacts.append(ExperimentalData())
                bioacts[-1].ligand_id = ligand_cache[row[' Ligand Name']].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = 'B'
                bioacts[-1].assay_description = None
                bioacts[-1].standard_activity_value  = round(float(row['ki Val']),2)
                bioacts[-1].p_activity_value = round(-math.log10(float(row['ki Val'])),2)
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = '='
                bioacts[-1].value_type = 'pKi'
                bioacts[-1].source = 'PDSP KiDatabase'
                bioacts[-1].document_chembl_id = None
                # BULK insert every X entries or last entry
                if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                    ExperimentalData.objects.bulk_create(bioacts)
                    print("Inserted", index, "out of", bio_entries, "bioactivities")
                    bioacts = []
