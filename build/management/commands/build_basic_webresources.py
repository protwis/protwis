from django.core.management.base import BaseCommand, CommandError
from common.models import WebResource
import logging

class Command(BaseCommand):
    
    logger = logging.getLogger(__name__)

    help = "Adding basic webresources: pubmed, pdb and iuphar"

    def handle(self, *args, **options):
        self.logger.info("ADDING BASIC WEB RESOURCES")
        #pubmed
        try:
            wr = WebResource.objects.get(slug='pubmed', name='NCBI PubMed')
        except WebResource.DoesNotExist:
            try:
                wr = WebResource.objects.create(slug='pubmed', name='NCBI PubMed', url='http://www.ncbi.nlm.nih.gov/pubmed/$index')
                self.logger.info("pubmed record successfully created")
            except Exception as msg:
                print(msg)
                self.logger.error("Couldn't create a record for pubmed")
        #pdb
        try:
            wr = WebResource.objects.get(slug='pdb', name='Protein Data Bank')
        except WebResource.DoesNotExist:
            try:
                wr = WebResource.objects.create(slug='pdb', name='Protein Data Bank', url='http://www.rcsb.org/pdb/explore/explore.do?structureId=$index')
                self.logger.info("pdb record successfully created")
            except Exception as msg:
                print(msg)
                self.logger.error("Couldn't create a record for pdb")
        #iuphar guidetopahrmacology
        try:
            wr = WebResource.objects.get(slug='iuphar', name='IUPHAR Guide to pharmacology')
        except WebResource.DoesNotExist:
            try:
                wr = WebResource.objects.create(slug='iuphar', name='IUPHAR Guide to pharmacology', url='http://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=$index&familyId=$family&familyType=GPCR')
                self.logger.info("iuphar record successfully created")
            except Exception as msg:
                print(msg)
                self.logger.error("Couldn't create a record for iuphar")
        self.logger.info("DONE")