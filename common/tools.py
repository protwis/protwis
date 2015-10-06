from django.conf import settings

import os
import yaml
import time
import logging
from urllib.request import urlopen
from urllib.error import HTTPError

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

def fetch_from_web_api(url):
    logger = logging.getLogger('build')
    logger.info('Fetching {}'.format(url))
    tries = 0
    max_tries = 5
    while tries < max_tries:
        if tries > 0:
            logger.warning('Failed fetching {}, retrying'.format(url))
        
        try:
            req = urlopen(url)
            return req
        except HTTPError:
            tries += 1
            time.sleep(2)
    
    logger.error('Failed fetching {} {} times, giving up'.format(url, max_tries))
    return False