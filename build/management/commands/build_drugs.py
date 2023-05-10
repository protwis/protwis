from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import IntegrityError

from common.models import WebResource, WebLink, Publication
from protein.models import Protein
from drugs.models import Drugs
from mutational_landscape.models import NHSPrescribings
from common.tools import test_model_updates

import pandas as pd
import os
import django.apps
import logging

class Command(BaseCommand):
    help = 'Build Drug and NHS Data'

    publication_cache = {}

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='store', dest='filename',
                            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

    # source file directory
    data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]

    def handle(self, *args, **options):
        try:
            self.purge_data()
            test_model_updates(self.all_models, self.tracker, initialize=True)
            self.create_drug_data()
            self.create_NHS()
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_data(self):
        try:
            Drugs.objects.all().delete()
            NHSPrescribings.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')
            self.logger.warning('Drugs mod not found: nothing to delete.')

    @staticmethod
    def fetch_publication(publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        if pd.isna(publication_doi) is True:
            return None

        if ("ISBN" in publication_doi) or (publication_doi == '0'):
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

        if publication_doi not in Command.publication_cache:
            try:
                wl = WebLink.objects.get(
                    index=publication_doi, web_resource__slug=pub_type)
            except WebLink.DoesNotExist:
                try:
                    wl = WebLink.objects.create(
                        index=publication_doi, web_resource=WebResource.objects.get(slug=pub_type))
                except IntegrityError:
                    wl = WebLink.objects.get(
                        index=publication_doi, web_resource__slug=pub_type)

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
                except Exception as e:
                    # if something off with publication, skip.
                    print("Build drugs Publication fetching error | module: fetch_publication. Row # is : " +
                          str(publication_doi) + ' ' + pub_type)
                    print(f'{type(e).__name__} {e}')

            Command.publication_cache[publication_doi] = pub
        else:
            pub = Command.publication_cache[publication_doi]

        return pub

    def read_csv_data(self, filename, file_suffix):

        if not filename:
            filename = next((fn for fn in os.listdir(self.data_dir) if fn.endswith(file_suffix)), None)

        filepath = os.sep.join([self.data_dir, filename])
        data = pd.read_csv(filepath, low_memory=False, encoding="ISO-8859-1")
        return data

    def create_drug_data(self, filename=False):
        print('CREATING DRUGDATA')

        data = self.read_csv_data(filename, 'drug_data.csv')
        data['PMID'] = data['PMID'].fillna('').astype(str)

        for _, row in data.iterrows():
            drugname = row['Drug Name'].split(",")[0]
            drugalias_raw = row['DrugAliases']
            drugalias = ['' if str(drugalias_raw) == 'nan' else str(drugalias_raw)][0]

            entry_name = row['EntryName']

            phase = row['Phase']
            PhaseDate = row['PhaseDate']
            ClinicalStatus = row['ClinicalStatus']
            moa = row['ModeOfAction']
            targetlevel = row['TargetCategory']

            drugtype = row['Drug Class']
            indication = row['Indication(s)'].title()
            novelty = row['Target_novelty']
            approval = row['Approval']
            status = row['Status']

            references = row['PMID']

            drug, _ = Drugs.objects.get_or_create(name=drugname,
                                                    synonym=drugalias,
                                                    drugtype=drugtype,
                                                    indication=indication,
                                                    novelty=novelty,
                                                    approval=approval,
                                                    phase=phase,
                                                    phasedate=PhaseDate,
                                                    clinicalstatus=ClinicalStatus,
                                                    moa=moa,
                                                    status=status,
                                                    targetlevel=targetlevel,
                                                    references=references)

            ref = references.split('|')
            try:
                for pmid in ref:
                    if pmid != '':
                        publication = Command.fetch_publication(pmid)
                        drug.publication.add(publication)
            except Exception as e:
                print(f'The Drugs {pmid} publication was not added to the data base')
                print(f'{type(e).__name__} {e} on build_drugs')

            # fetch protein
            try:
                p = Protein.objects.get(entry_name=entry_name)
                drug.target.add(p)
            except Protein.DoesNotExist:
                print('Protein not found for entry_name {}'.format(entry_name))

            drug.save()

        self.logger.info('COMPLETED CREATING DRUGDATA')

    def create_NHS(self, filename=False):
        print('Creating NHS')

        nhs_data = self.read_csv_data(filename, 'nhs.csv')

        for _, entry in nhs_data.iterrows():
            date = entry['date']
            quantity = entry['quantity']
            items = entry['items']
            actual_cost = entry['actual_cost']
            drugCode = entry['drugCode']
            op_name = entry['drugName']
            bnf_section_raw = entry['section']
            bnf_section_name = bnf_section_raw.split(': ')[1]
            drugNameQuery = entry['drugNameQuery']
            try:
                drugname = Drugs.objects.filter(name=drugNameQuery)[0]
            except Exception:
                self.logger.warning('Drug not found for chemical {}'.format(drugNameQuery))
                continue

            NHSPrescribings.objects.get_or_create(date=date,
                                                    quantity=quantity,
                                                    items=items,
                                                    actual_cost=actual_cost,
                                                    drugCode=drugCode,
                                                    op_name=op_name,
                                                    bnf_section=bnf_section_name,
                                                    drugname=drugname)

        self.logger.info('COMPLETED CREATING NHS PRESCRIBINGS')
