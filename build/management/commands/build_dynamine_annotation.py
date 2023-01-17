from django.core.cache import cache
from django.core.management.base import BaseCommand
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify

from django.conf import settings
from common.tools import test_model_updates
from residue.models import Residue, ResidueDataType, ResidueDataPoint
from protein.models import Protein, ProteinSet

import logging
from urllib import request, parse
import json
import time
import os
import django.apps

class Command(BaseBuild):
    help = "Add dynamine annotations."

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        parser.add_argument("-p", "--proc",
                            type=int,
                            action="store",
                            dest="proc",
                            default=1,
                            help="Number of processes to run")

    def handle(self, *args, **options):
        # Purge previous dynamine data
        ResidueDataPoint.objects.filter(data_type__slug=slugify("dynamine")).delete()
        ResidueDataType.objects.filter(slug=slugify("dynamine"), name="Dynamine Prediction").delete()

        # Create a new dynamine residue data type here (previously a parallel unsafe operation)
        dynamine_type = ResidueDataType(slug=slugify("dynamine"), name="Dynamine Prediction")
        dynamine_type.save()

        # All human proteins and xtaled
        self.proteins = list(set(list(Protein.objects.filter(sequence_type__slug='wt',species__common_name="Human").all())+list(ProteinSet.objects.get(name="All").proteins.all())))
        #self.proteins = list( set(list(ProteinSet.objects.get(name="All").proteins.all())))
        self.prepare_input(options["proc"], self.proteins)
        # self.logger.info('Finishing adding dynamine annotations')
        test_model_updates(self.all_models, self.tracker, check=True)

    def get_dynamine_prediction(self, protein):
        # DEPRECATED First: try the local Django cache and return if exists
        # dynamine_results = cache.get("dynamine_prediction_%s" % (protein.entry_name))
        # if dynamine_results:
        #     return dynamine_results

        # Second: try the dynamine file cache and return if exists
        dynamine_dir = os.sep.join([settings.DATA_DIR, 'structure_data','dynamine_cache'])
        cache_file = os.sep.join([dynamine_dir, "{}_cache.json".format(protein.entry_name)])
        if os.path.isfile(cache_file):
            with open(cache_file, 'r') as openfile:
                # Grab from filecache
                dynamine_results = json.load(openfile)

                # Write to DJANGO cache
                cache.set('dynamine_prediction_%s' %
                          (protein.entry_name), dynamine_results, 60 * 60 * 24 * 7)

                # Return results
                return dynamine_results

        # Third: request from the dynamine server
        json_api_key = "26cf5c434a171cab9220666030cd981bdbe485a449729a7c8c6272b9"
        job = {"tool_list": ["dynamine"],
               "TOKEN": json_api_key,
               protein.entry_name: protein.sequence,
               }
        url = 'https://bio2byte.be/msatools/api/'
        response_code = 0
        try:
            # "POST" the job information to the DynaMine API
            data = parse.urlencode(job).encode()
            req = request.Request(url, data=data)
            reqr = request.urlopen(req)
            resp = json.loads(reqr.read().decode("UTF-8"))
        except Exception as e:
            resp = {}
            print("Error starting job", e)

        hash_id = resp["hash_id"]
        tries = 0
        response_code = -1
        while response_code != "200":
            time.sleep(5)
            tries += 1
            try:
                req = request.Request("https://bio2byte.be/msatools/api/queue/" + hash_id)
                resp = json.loads(request.urlopen(req).read().decode('UTF-8'))
                if "status" in resp and resp["status"] not in [200, 202]:
                    print("Errors when querying DynaMine server")
                    exit(0)
                elif "results" in resp:
                    response_code = "200"
                    # Error in new DynaMine API => tries processing tool_list and TOKEN entries
                    resp["results"] = {"predictions": { resp["results"][-1]["proteinID"] : resp["results"][-1]["backbone"]}}

            except Exception as e:
                resp = {}
                print("Error polling job", e)
            if tries > 60:
                break

        # Write new results to DJANGO cache
        cache.set('dynamine_prediction_%s' %
                  (protein.entry_name), resp, 60 * 60 * 24 * 7)  # 7 days

        # Write new results to file cache
        with open(cache_file, 'w+') as f:
            f.write(json.dumps(resp))
            f.close()

        return resp

    def save_dynamine_prediction(self, protein):
        r = self.get_dynamine_prediction(protein)
        dynamine, created = ResidueDataType.objects.get_or_create(
            slug=slugify("dynamine"), name="Dynamine Prediction")
        residues = Residue.objects.filter(
            protein_conformation__protein=protein).all()

        c = 0
        if "results" in r and protein.entry_name in r["results"]["predictions"]:
            predictions = r["results"]["predictions"][protein.entry_name]
            bulk = []
            for i, p in enumerate(predictions):
                # fetch residue
                try:
                    r = residues.filter(sequence_number=i + 1).get()
                    if isinstance(p, float):
                        data = p # process results > Feb-2022
                    else:
                        data = p[1] # process results <= Feb-2022

                    c += 1
                    bulk.append(ResidueDataPoint(data_type=dynamine, residue=r, value=data))

                except:
                    print("Missing residue for", protein.entry_name, i + 1)
            ResidueDataPoint.objects.bulk_create(bulk)
        else:
            print(protein, r)
            print("ERROR processing DynaMine for", protein.entry_name)
            exit(0)

        #print("PROCESSED", protein, "residue count", c, "vs", len(protein.sequence))
        return c

    # @transaction.atomic
    def main_func(self, positions, iteration, count, lock):
        # Process all proteins
        while count.value < len(self.proteins):
            with lock:
                p = self.proteins[count.value]
                count.value += 1
                self.logger.info("Generating dynamine data for '{}'... ({} out of {})".format(
                    p, count.value, len(self.proteins)))
            # if 'opsd_bovin'!=str(p):
            #     continue
            dynamine = self.save_dynamine_prediction(p)
