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
                wr = WebResources.objects.create(slug='pubmed', name='NCBI PubMed', link='http://www.ncbi.nlm.nih.gov/pubmed/$index')
                self.logger.info("pubmed record succesfully created")
            except msg:
                print(msg)
                self.logger.error("Couldn't create a record for pubmed")
        #pdb
        try:
            wr = WebResource.objects.get(slug='pdb', name='Protein Data Bank')
        except WebResource.DoesNotExist:
            try:
                wr = WebResources.objects.create(slug='pdb', name='Protein Data Bank', link='http://www.rcsb.org/pdb/explore/explore.do?structureId=$index')
                self.logger.info("pdb record succesfully created")
            except msg:
                print(msg)
                self.logger.error("Couldn't create a record for pubmed")
        #iuphar guidetopahrmacology
        try:
            wr = WebResource.objects.get(slug='iuphar', name='IUPHAR Guide to pharmacology')
        except WebResource.DoesNotExist:
            try:
                wr = WebResources.objects.create(slug='iuphar', name='IUPHAR Guide to pharmacology', link='http://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=$index&familyId=$family&familyType=GPCR')
                self.logger.info("iuphar record succesfully created")
            except msg:
                print(msg)
                self.logger.error("Couldn't create a record for pubmed")
        