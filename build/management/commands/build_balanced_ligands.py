from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein
from ligand.models import  Ligand, BiasedData, BalancedLigands

from common.models import Publication
import logging
import math
import pandas as pd
import numpy as np

#Hello and welcome to the build of the Balanced Ligands.
#This build is totally subordinated to the data in the BiasedData Model
#and exploits an adjusted version of the OnTheFly function to parse
#the data and calculate those ligands that show a balanced preference
#in terms of selection between two (or more) activation pathways.
#These data are not manually curated but are artifically generated
#through calculation based on the assumption that, if the ratio of
#Log(Tau/KA) or Log(Emax/EC50) between two pathways is similar to 1
#to a degree of extent of 20%, they show same prefere towards both pathways,
#though can be defined as 'balanced' between those pathways.

class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    publication_cache = {}
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('BalancedLigands.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)

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
                            help='Purge existing balanced ligands records')
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run',
                            default=False)

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        if options['purge']:
            try:
                Command.purge_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        Command.model_assemble()

    @staticmethod
    def purge_data():
        print("\n# Purging data")
        balanced_data = BalancedLigands.objects.all()
        balanced_data.delete()
        print("\n# Old data removed")

    @staticmethod
    def model_assemble():
        """
        Fetch data to models
        Saves to DB
        """
        print('*** Starting *** \n')
        print('\n#1 Parsing data from BiasedData model')
        receptors = list(BiasedData.objects.all().values_list("receptor_id", flat=True).distinct())
        queried_db = Command.parse_data(receptors)
        print('\n#2 Generating the reference balanced ligands dataframe')
        balanced_db = Command.creating_balanced_database(queried_db)
        print('\n#3 Pushing the generated dataframe to the model')
        Command.create_model(balanced_db)
        print('\n *** Finished! ***')

    @staticmethod
    def define_ligand_pathways(master):
        ligands = {}
        for key in master:
            if master[key]['ligand_id'] not in ligands.keys():
                ligands[master[key]['ligand_id']] = []
            ligands[master[key]['ligand_id']].append(key)
        return ligands

    @staticmethod
    def assess_pathway_preferences(comparisons, tested, subtype=False):
        pathway_preference = {}
        #calculate values for ranking (or replace with qualitative activity when missing)
        for assay in comparisons.keys():
            for test in comparisons[assay]:
                if tested[test]['ligand_id'] not in pathway_preference.keys():
                    pathway_preference[tested[test]['ligand_id']] = {}
                if subtype:
                    path_label = str(tested[test]['primary_effector_family']) + ' - ' + str(tested[test]['primary_effector_subtype'])
                else:
                    path_label = tested[test]['primary_effector_family']
                Tau_KA = tested[test]['Tau_KA']
                try:
                    Emax_EC50 = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                except (TypeError, ValueError):
                    Emax_EC50 = None
                pathway_preference[tested[test]['ligand_id']][path_label] = [Tau_KA, Emax_EC50]
        #ranking accordingly to Log(Tau/KA) or Log(Emax/EC50) (depending on how many missing values are)
        for val in pathway_preference.keys():
            temp = []
            none_tau = len([pathway_preference[val][key][0] for key in pathway_preference[val] if pathway_preference[val][key][0] is None])
            none_emax = len([pathway_preference[val][key][1] for key in pathway_preference[val] if pathway_preference[val][key][1] is None])
            if none_tau <= none_emax:
                pathway_preference[val] = list(dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if item[1][0] == None else item[1][0], reverse=True)).keys())
            else:
                pathway_preference[val] = list(dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if item[1][1] == None else item[1][1], reverse=True)).keys())
        #provide ranked keys
            for path in pathway_preference[val]:
                if subtype:
                    if path.split(' - ')[1] == 'None':
                        temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path.split(' - ')[0] and tested[key]['primary_effector_subtype'] ==  None)][0])
                    else:
                        temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path.split(' - ')[0] and tested[key]['primary_effector_subtype'] == path.split(' - ')[1])][0])
                else:
                    temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path)][0])
            pathway_preference[val] = temp
        return pathway_preference

    @staticmethod
    def calculate_delta(comparisons, tested, subtype=False):
        ranking = Command.assess_pathway_preferences(comparisons, tested, subtype)
        #STEPS
        calculated_output = {}
        for drug in ranking.keys():
            #perform analysis only when we have multiple pathways to actually compare
            if len(ranking[drug]) > 1:
            #Set reference pathway (first in list)
                path_count = 0
                while path_count < len(ranking[drug]) - 1:
                    reference = ranking[drug][path_count]
                    for test in ranking[drug][path_count+1:]:
                        path_count +=1
                        #Match pathway levels + skip matching if Arrestin involved
                        if (tested[test]['pathway_level'] == tested[reference]['pathway_level']) or ('Arrestin' in [tested[reference]['primary_effector_family'], tested[test]['primary_effector_family']]):
                            if subtype:
                                #make a check for excluding missing subtype info,
                                #because that way we are not investigating subtype properly
                                if (tested[reference]['primary_effector_subtype'] == None) or (tested[test]['primary_effector_subtype'] == None):
                                    continue
                                Pathway_Pair = tested[reference]['primary_effector_subtype'] + ' - ' + tested[test]['primary_effector_subtype']
                            else:
                                Pathway_Pair = tested[reference]['primary_effector_family'] + ' - ' + tested[test]['primary_effector_family']

                            ID = str(drug) + '_' + Pathway_Pair
                            calculated_output[ID] = {}
                            calculated_output[ID]['receptor_id'] = tested[test]['receptor_id']
                            calculated_output[ID]['ligand_id'] = tested[test]['ligand_id']
                            calculated_output[ID]['doi'] = tested[test]['doi']
                            calculated_output[ID]['Comparison'] = Pathway_Pair
                            try:
                                delta_logtauka = round(tested[reference]['Tau_KA'] - tested[test]['Tau_KA'], 3)
                            except TypeError:
                                delta_logtauka = None
                            try:
                                tested[reference]['Log(Emax/EC50)'] = round(math.log((tested[reference]['Emax']/tested[reference]['EC50']),10), 3)
                            except (TypeError, ValueError):
                                tested[reference]['Log(Emax/EC50)'] = None
                            try:
                                tested[test]['Log(Emax/EC50)'] = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                            except (TypeError, ValueError):
                                tested[test]['Log(Emax/EC50)'] = None
                            try:
                                delta_logemaxec50 = round(tested[reference]['Log(Emax/EC50)'] - tested[test]['Log(Emax/EC50)'], 3)
                            except TypeError:
                                delta_logemaxec50 = None

                            calculated_output[ID]['Delta Log(Tau/KA)'] = delta_logtauka
                            calculated_output[ID]['Delta Log(Emax/EC50)'] = delta_logemaxec50
                            if subtype:
                                calculated_output[ID]['subtype'] = 'YES'
                            else:
                                calculated_output[ID]['subtype'] = 'NO'

        return calculated_output

    @staticmethod
    def AdjustedFly(receptor_id, subtype=False):
        #fetching data given the receptor id
        test_data = BiasedData.objects.filter(receptor=receptor_id, Emax__gte=90)
        pub_ids = list(BiasedData.objects.filter(receptor=receptor_id).values_list("publication", flat=True).distinct())

        # Performance: first collect all publication and ligand data
        pub_objs = Publication.objects.filter(id__in=pub_ids).values_list("id", "web_link_id__index", "year", "journal_id__name", "authors")
        pub_objs_dict = {pub_obj[0]:pub_obj[1:] for pub_obj in pub_objs}

        publications = {}
        for entry in test_data:
            if entry.publication_id not in publications.keys():
                publications[entry.publication_id] = {}
            if entry.id not in publications[entry.publication_id].keys():
                if entry.publication_id in pub_objs_dict:
                    pub_data = pub_objs_dict[entry.publication_id]
                else:
                    pub_data = Publication.objects.filter(id=entry.publication_id).values_list("web_link_id__index")

                publications[entry.publication_id][entry.id] = {'doi': pub_data[0],
                                                                'receptor_id': entry.receptor_id,
                                                                'ligand_id': entry.ligand_id,
                                                                'primary_effector_family': entry.primary_effector_family,
                                                                'primary_effector_subtype': entry.primary_effector_subtype,
                                                                'pathway_level': entry.pathway_level,
                                                                'EC50': entry.EC50,
                                                                'Emax': entry.Emax,
                                                                'Tau_KA': entry.Tau_KA,}

        for pub in list(publications.keys()):
            ligands = Command.define_ligand_pathways(publications[pub])
            publications[pub] = Command.calculate_delta(ligands, publications[pub], subtype)

        return publications

    @staticmethod
    def parse_data(receptors):
        pathway_dump = pd.DataFrame()
        for protein in receptors:
            data = Command.AdjustedFly(protein)
            data_subtype = Command.AdjustedFly(protein, True)
            for publication in data.keys():
                for row in data[publication]:
                    pathway_dump = pathway_dump.append(data[publication][row], ignore_index=True)
            for publication in data_subtype.keys():
                for row in data_subtype[publication]:
                    pathway_dump = pathway_dump.append(data_subtype[publication][row], ignore_index=True)
        return pathway_dump

    @staticmethod
    def creating_balanced_database(parsed_data):
        parsed_data['preferred_metric'] = parsed_data['Delta Log(Tau/KA)'].combine_first(parsed_data['Delta Log(Emax/EC50)'])
        filtered = parsed_data.loc[~parsed_data['preferred_metric'].isnull()]

        test = filtered[filtered['preferred_metric'].between(-0.2, 0.2, inclusive="neither")]

        balanced_db = test[['receptor_id', 'ligand_id', 'doi', 'Comparison', 'Delta Log(Tau/KA)', 'Delta Log(Emax/EC50)', 'preferred_metric', 'subtype']].copy()
        balanced_db[['Pathway 1', 'Pathway 2']] = balanced_db['Comparison'].str.split(' - ', 1, expand=True)
        #this step is needed to undiscriminate order of pathways in comparisons
        #while assessing the balance reference ligand for that specific pathway pair
        balanced_db['sorter'] = ''
        for index, row in balanced_db.iterrows():
            row['sorter'] = ''.join(list(set([row['Pathway 1'], row['Pathway 2']])))

        recs = list(balanced_db['receptor_id'].unique())
        balanced_db['Ranking'] = np.nan
        to_be_parsed = balanced_db.groupby(['receptor_id', 'sorter', 'doi']).size().reset_index().rename(columns={0:'count'})
        the_tuples = list(zip(to_be_parsed.receptor_id, to_be_parsed.sorter, to_be_parsed.doi))

        for pair in the_tuples:
            piece = balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2])]
            if len(piece) != 1:
                ranked_dict = dict(zip(piece.ligand_id, abs(piece.preferred_metric)))
                ranked_dict = {ligand: diff_value for ligand, diff_value in sorted(ranked_dict.items(), key=lambda item: item[1])}
                Rank = 1
                for key in ranked_dict.keys():
                    balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2]) & (balanced_db['ligand_id'] == key), 'Ranking'] =  Rank
                    Rank += 1
            else:
                balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2]), 'Ranking'] =  1

        return balanced_db

    @staticmethod
    def create_model(balanced_db):
        reference_balanced_db = balanced_db.loc[balanced_db['Ranking'] == 1].drop_duplicates()
        for index, row in reference_balanced_db.iterrows():
            p = Command.fetch_protein(row['receptor_id'])
            l = Command.fetch_ligand(row['ligand_id'])
            pub = Command.fetch_publication(row['doi'])
            subtype = True if row['subtype'] == 'YES' else False
            balanced_data = BalancedLigands(
                        ligand = l,   #link to ligand model
                        receptor = p, #link to protein model
                        first_pathway = row['Pathway 1'],
                        second_pathway = row['Pathway 2'],
                        delta_logEmaxEC50 = row['Delta Log(Emax/EC50)'],
                        delta_logTauKA = row['Delta Log(Tau/KA)'],
                        publication = pub,
                        subtype_balanced = subtype,
                        )
            balanced_data.save()


    @staticmethod
    def fetch_protein(protein_id):
        """
        fetch receptor with Protein model
        requires: protein id.
        Should NEVER return None given the
        logic of this build script
        """
        try:
            test = None
            if Protein.objects.filter(id=protein_id):
                protein = Protein.objects.filter(id=protein_id)
                test = protein.get()
            return test
        except:
            return None

    @staticmethod
    def fetch_ligand(ligand_id):
        """
        fetch ligands with Ligand model
        requires: ligand id.
        Should NEVER return None given the
        logic of this build script
        """
        try:
            test = None
            if Ligand.objects.filter(id=ligand_id):
                drug = Ligand.objects.filter(id=ligand_id)
                test = drug.get()
            return test
        except:
            return None

    @staticmethod
    def fetch_publication(publication_doi):
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
        try:
            if publication_doi not in Command.publication_cache:
                pub = False
                if pub_type == 'doi':
                    pub = Publication.get_or_create_from_doi(publication_doi)
                elif pub_type == 'pubmed':
                    pub = Publication.get_or_create_from_pubmed(publication_doi)
                Command.publication_cache[publication_doi] = pub
            else:
                pub = Command.publication_cache[publication_doi]
        except:
            pub = Publication.objects.filter(web_link__index = publication_doi).first()
        return pub
