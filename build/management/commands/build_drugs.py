from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from common.models import WebResource, WebLink, Publication
from protein.models import Protein
from drugs.models import Drugs


from optparse import make_option
import logging
import csv
import os
import pandas as pd

class Command(BaseCommand):
    help = 'Build Drug Data'
    publication_cache = {}

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

        # source file directory
    drugdata_data_dir = os.sep.join([settings.DATA_DIR, 'drug_data']) #settings.DATA_DIR

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.purge_drugs()
            self.create_drug_data()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_drugs(self):
        try:
            Drugs.objects.all().delete()
        except Drugs.DoesNotExist:
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
                except:
                    # if something off with publication, skip.
                    print("Publication fetching error | module: fetch_publication. Row # is : " +
                          str(publication_doi) + ' ' + pub_type)

            Command.publication_cache[publication_doi] = pub
        else:
            pub = Command.publication_cache[publication_doi]

        return pub
    

    def create_drug_data(self, filenames=False):
        print('CREATING DRUGDATA')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.drugdata_data_dir) if fn.endswith('drug_data.csv')]

        for filename in filenames:

            filepath = os.sep.join([self.drugdata_data_dir, filename])

            data = pd.read_csv(filepath, low_memory=False, encoding = "ISO-8859-1", dtype={'PMID': str})
            data['PMID'] = data['PMID'].fillna('')

            for _, row in data.iterrows():
                drugname = row['Drug Name'].split(",")[0]
                trialname = row['Trial name']
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

                # fetch protein

                drug, created = Drugs.objects.get_or_create(name=drugname, 
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
                                                            )
                try:
                    p = Protein.objects.get(entry_name=entry_name)
                    drug.target.add(p)
                except Protein.DoesNotExist:
                    print('Protein not found for entry_name {}'.format(entry_name))
                    continue

                ref = references.split('|')

                try:
                    for pmid in ref:
                        publication = Command.fetch_publication(pmid)
                        drug.publication.add(publication)
                except:
                    publication = None


                drug.save()

                # target_list = drug.target.all()

        self.logger.info('COMPLETED CREATING DRUGDATA')
