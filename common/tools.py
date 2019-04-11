from django.conf import settings
from django.utils.text import slugify
from django.core.cache import cache

import os
import yaml
import time
import logging
import urllib
from urllib.parse import quote
from urllib.request import urlopen
from urllib.error import HTTPError
import json
import gzip
from io import BytesIO
from string import Template
from Bio import Entrez, Medline
import xml.etree.ElementTree as etree 



def save_to_cache(path, file_id, data):
    create_cache_dirs(path)
    cache_dir_path = os.sep.join([settings.BUILD_CACHE_DIR] + path)
    cache_file_path = os.sep.join([cache_dir_path, file_id + '.yaml'])
    with open(cache_file_path, 'w') as cache_file:
        yaml.dump(data, cache_file, default_flow_style=False)

def fetch_from_cache(path, file_id):
    cache_dir_path = os.sep.join([settings.BUILD_CACHE_DIR] + path)
    cache_file_path = os.sep.join([cache_dir_path, file_id + '.yaml'])
    if os.path.isfile(cache_file_path):
        with open(cache_file_path) as cache_file:
            return yaml.load(cache_file)
    else:
        return False

def create_cache_dirs(path):
    cache_file_path = os.sep.join([settings.BUILD_CACHE_DIR] + path)
    os.makedirs(cache_file_path, exist_ok=True)
    intermediate_path = settings.BUILD_CACHE_DIR
    for directory in path:
        intermediate_path = os.sep.join([intermediate_path, directory])
        os.chmod(intermediate_path, 0o777)

def fetch_from_web_api(url, index, cache_dir=False, xml=False, raw=False):
    logger = logging.getLogger('build')

    # slugify the index for the cache filename (some indices have symbols not allowed in file names (e.g. /))
    index_slug= slugify(index)
    cache_file_path = '{}/{}'.format('/'.join(cache_dir), index_slug)
    
    # try fetching from cache
    if cache_dir:
        d = cache.get(cache_file_path)
        # d = fetch_from_cache(cache_dir, index_slug)
        if not d is None:
            logger.info('Fetched {} from cache'.format(cache_file_path))
            return d
    
    # if nothing is found in the cache, use the web API
    full_url = Template(url).substitute(index=quote(str(index), safe=''))
    logger.info('Fetching {}'.format(full_url))
    tries = 0
    max_tries = 5
    while tries < max_tries:
        if tries > 0:
            logger.warning('Failed fetching {}, retrying'.format(full_url))
        
        try:
            req = urlopen(full_url)
            if full_url[-2:]=='gz' and xml:
                try:
                    buf = BytesIO( req.read())
                    f = gzip.GzipFile(fileobj=buf)
                    data = f.read()
                    d = etree.fromstring(data)
                except:
                    return False
            elif xml:
                try:
                    d = etree.fromstring(req.read().decode('UTF-8'))
                except:
                    return False
            elif raw:
                try:
                    d = req.read().decode('UTF-8')
                except:
                    return False
            else:
                d = json.loads(req.read().decode('UTF-8'))
        except HTTPError as e:
            tries += 1
            if e.code == 404:
                logger.warning('Failed fetching {}, 404 - does not exist'.format(full_url))
                tries = max_tries #skip more tries
                return False
            elif e.code == 400:
                logger.warning('Failed fetching {}, 400 - does not exist'.format(full_url))
                tries = max_tries #skip more tries
                return False
            else:
                time.sleep(2)
        except urllib.error.URLError as e:
            # Catches 101 network is unreachable -- I think it's auto limiting feature
            tries +=1 
            time.sleep(2)
        else:
            # save to cache
            if cache_dir:
                # save_to_cache(cache_dir, index_slug, d)
                cache.set(cache_file_path, d, 60*60*24*7) #7 days
                logger.info('Saved entry for {} in cache'.format(cache_file_path))
            return d
    
    # give up if the lookup fails 5 times
    logger.error('Failed fetching {} {} times, giving up'.format(full_url, max_tries))
    return False

def fetch_from_entrez(index, cache_dir=False):
    logger = logging.getLogger('build')

    # slugify the index for the cache filename (some indices have symbols not allowed in file names (e.g. /))
    index_slug= slugify(index)
    cache_file_path = '{}/{}'.format('/'.join(cache_dir), index_slug)

    # try fetching from cache
    if cache_dir:
        d = fetch_from_cache(cache_dir, index_slug)
        if d:
            logger.info('Fetched {} from cache'.format(cache_file_path))
            return d
    
    # if nothing is found in the cache, use the web API
    logger.info('Fetching {} from Entrez'.format(index))
    tries = 0
    max_tries = 5
    while tries < max_tries:
        if tries > 0:
            logger.warning('Failed fetching pubmed {}, retrying'.format(str(index)))
            
        try:
            Entrez.email = 'info@gpcrdb.org'
            handle = Entrez.efetch(
                db="pubmed", 
                id=str(index), 
                rettype="medline", 
                retmode="text"
            )
        except:
            tries += 1
            time.sleep(2)
        else:
            d = Medline.read(handle)

            # save to cache
            save_to_cache(cache_dir, index_slug, d)
            logger.info('Saved entry for {} in cache'.format(cache_file_path))
            return d