from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import IntegrityError

from common.models import WebResource, WebLink, Publication
from protein.models import TissueExpression
from protein.models import Protein
from mutational_landscape.models import NHSPrescribings
from common.tools import test_model_updates

import pandas as pd
import os
import django.apps
import logging





class Command(BaseCommand):
    help = 'Build Tissue Expression Data'

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
            TissueExpression.objects.all().delete()
            #NHSPrescribings.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')
            self.logger.warning('TissueExpression mod not found: nothing to delete.')

    @staticmethod
    def read_csv_data(filename):
        filepath = os.sep.join([data_dir, filename])
        data = pd.read_csv(filepath, low_memory=False,encoding="ISO-8859-1")
        return data

    # for testing #
    
    Filename = '08_TargetPrioritazion_AllData.csv'
    Filename2 = "03_FINAL_DATA_UPDATED.csv"
    data = read_csv_data(Filename)
    data2 = read_csv_data(Filename)
    # Cleaning the tissue data #

    Tissue_cols = ['entry_name'] + [col for col in data if col.startswith('Tissue')]
    
    data = data[Tissue_cols].drop_duplicates().dropna(subset=Tissue_cols[1:])
    
    # Function for renaming the columns #
    def renaming_cols(col):
        if col.startswith('Tissue'):
            keep_name = col.split(' - ')[1].split(' [')[0].replace(' 1','').replace(' ','_')
            return keep_name
        else:
            pass
    # Rename the columns #
    data.rename(columns=renaming_cols,inplace=True)
    
    ## Push and create data for tissue expression ##

    def fetch_protein(target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            protein = Protein.objects.get(entry_name=target)
            return protein
        except:
            print('No protein found for this entry name')
            return None

    def create_Expression_data(self, filename=False):
        print('CREATING EXPRESSION DATA')

        #data = self.read_csv_data(filename, '08_TargetPrioritazion_AllData.csv')

    for i, row in data.iterrows():
        drug, i = TissueExpression.objects.get_or_create(
        protein = fetch_protein(row['entry_name']),
        adipose_tissue = row['adipose_tissue'],
        adrenal_gland = row['adrenal_gland'],
        amygdala = row['amygdala'],
        appendix = row['appendix'],
        basal_ganglia = row['basal_ganglia'],
        bone_marrow = row['bone_marrow'],
        breast = row['breast'],
        cerebellum = row['cerebellum'],
        cerebral_cortex = row['cerebral_cortex'],
        cervix = row['cervix'],
        choroid_plexus = row['choroid_plexus'],
        colon = row['colon'],
        duodenum = row['duodenum'],
        endometrium = row['endometrium'],
        epididymis = row['epididymis'],
        esophagus = row['esophagus'],
        fallopian_tube = row['fallopian_tube'],
        gallbladder = row['gallbladder'],
        heart_muscle = row['heart_muscle'],
        hippocampal_formation = row['hippocampal_formation'],
        hypothalamus = row['hypothalamus'],
        kidney = row['kidney'],
        liver = row['liver'],
        lung = row['lung'],
        lymph_node = row['lymph_node'],
        midbrain = row['midbrain'],
        ovary = row['ovary'],
        pancreas = row['pancreas'],
        parathyroid_gland = row['parathyroid_gland'],
        pituitary_gland = row['pituitary_gland'],
        placenta = row['placenta'],
        prostate = row['prostate'],
        rectum = row['rectum'],
        retina = row['retina'],
        salivary_gland = row['salivary_gland'],
        seminal_vesicle = row['seminal_vesicle'],
        skeletal_muscle = row['skeletal_muscle'],
        skin = row['skin'],
        small_intestine = row['small_intestine'],
        smooth_muscle = row['smooth_muscle'],
        spinal_cord = row['spinal_cord'],
        spleen = row['spleen'],
        stomach = row['stomach'],
        testis = row['testis'],
        thymus = row['thymus'],
        thyroid_gland = row['thyroid_gland'],
        tongue = row['tongue'],
        tonsil = row['tonsil'],
        urinary_bladder = row['urinary_bladder'],
        vagina = row['vagina'])
        drug.save()

        self.logger.info('COMPLETED CREATING EXPRESSION DATA')


## console script: ##
        
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import IntegrityError

from common.models import WebResource, WebLink, Publication
from protein.models import TissueExpression
from protein.models import Protein
from mutational_landscape.models import NHSPrescribings
from common.tools import test_model_updates

import pandas as pd
import os
import django.apps
import logging

data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])

def read_csv_data(filename):
filepath = os.sep.join([data_dir, filename])
data = pd.read_csv(filepath, low_memory=False,encoding="ISO-8859-1")
return data

Filename = '08_TargetPrioritazion_AllData.csv'

data = read_csv_data(Filename)