from django.core.management.base import BaseCommand
from build.management.commands.base_build import Command as BaseBuild

from residue.models import ResidueDataType, ResidueDataPoint
from protein.models import *


import logging
from urllib import request, parse
import json,time

class Command(BaseBuild):
    help = 'Add dynamine annotations.'

    logger = logging.getLogger(__name__)
    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')

    def handle(self, *args, **options):

        # All human proteins and xtaled
        # self.proteins = list(set(list(Protein.objects.filter(sequence_type__slug='wt',species__common_name="Human").all())+list(ProteinSet.objects.get(name='All').proteins.all())))
        self.proteins = list(set(list(ProteinSet.objects.get(name='All').proteins.all())))
        # print(self.proteins)

        self.prepare_input(options['proc'], self.proteins)
        # self.logger.info('Finishing adding dynamine annotations')
    def get_dynamine_prediction(self, protein):
        json_api_key = '26cf5c434a171cab9220666030cd981bdbe485a449729a7c8c6272b9'
        job = {'protocol': '1.0',
               'json_api_key': json_api_key,
               'sequences': {protein.entry_name: protein.sequence},
               'predictions_only': 'true',
        }
        url = 'http://dynamine.ibsquare.be/batch_request'

        d = cache.get('dynamine_prediction_%s' % (protein.entry_name))
        if d:
            return d
        try:
            data = parse.urlencode({'batch':json.dumps(job)}).encode()
            req =  request.Request(url, data=data) # this will make the method "POST"
            resp = json.loads(request.urlopen(req).read().decode('UTF-8'))
        except Exception as e:
            resp = {}
            print('error starting job',e)
        job_id = resp['job_id']
        poll = {'protocol': '1.0',
               'json_api_key': json_api_key,
               'job_id': job_id
        }
        resp['status'] = 'Started'
        tries=0
        while resp['status'] != 'completed':
            time.sleep(2)
            tries += 1
            try:
                data = parse.urlencode({'batch':json.dumps(poll)}).encode()
                req =  request.Request(url, data=data) # this will make the method "POST"
                resp = json.loads(request.urlopen(req).read().decode('UTF-8'))
            except Exception as e:
                resp = {}
                print('error polling job',e)
            if tries>60:
                break

        cache.set('dynamine_prediction_%s' % (protein.entry_name), resp, 60*60*24*7) #7 days
        return resp

    def save_dynamine_prediction(self,protein):
        r = self.get_dynamine_prediction(protein)
        dynamine, created = ResidueDataType.objects.get_or_create(slug=slugify('dynamine'), name='Dynamine Prediction')
        residues = Residue.objects.filter(protein_conformation__protein=protein).all()
        c = r['status']
        if r['status'] == 'completed':
            predictions = r['results']['predictions'][protein.entry_name]
            # print(predictions)
            c = 0
            for i, p in enumerate(predictions):
                # fetch residue
                try:
                    r = residues.filter(sequence_number=i+1).get()
                    # print(protein,r,p[1],r.pk,i)
                    point, created = ResidueDataPoint.objects.get_or_create(data_type=dynamine, residue=r, value=p[1])
                    if created:
                        c += 1
                except:
                    print("Missing residue for",protein.entry_name,i+1)
        return c


    # @transaction.atomic
    def main_func(self, positions, iteration,count,lock):
        while count.value<len(self.proteins):
            with lock:
                p = self.proteins[count.value]
                count.value +=1 
                self.logger.info('Generating dynamine data for \'{}\'... ({} out of {})'.format(p, count.value, len(self.proteins)))
            # if 'opsd_bovin'!=str(p):
            #     continue
            dynamine = self.save_dynamine_prediction(p)
            # print(p,dynamine)

