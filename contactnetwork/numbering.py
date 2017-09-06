import json
import urllib.request
import os


def get_sequence_to_gpcrdb_map(pdb_name):
    # Download the map as JSON
    response = urllib.request.urlopen('http://test.gpcrdb.org/services/residues/%s/?format=json' % pdb_name)

    # Check that the request was succesful
    if response.getcode() != 200:
        raise Exception('Could not retrieve residue numbering scheme from GPCRdb.')

    # Check if there are any entries
    if not response:
        raise Exception('No numbering scheme for PDB file %s in GPCRdb.' % pdb_name)

    # Get the contents
    contents = response.read()

    # Parse the response as JSON
    json_list = json.loads(contents.decode("utf-8"))

    # The map
    seq_to_gpcr_map = {e[u'sequence_number']: e[u'display_generic_number'] for e in json_list}

    seq_to_gpcr_gen_map = {}

    # Convert into generic position, i.e. 3.28x29 -> 3x29
    for seq, gen in seq_to_gpcr_map.items():
        if gen is not None:
            tokens = gen.split('.')
            tm = tokens[0]
            tokens = tokens[1].split('x')
            gen = tokens[1]
            seq_to_gpcr_gen_map[seq] = tm + 'x' + gen
        else:
            seq_to_gpcr_gen_map[seq] = gen

    return seq_to_gpcr_gen_map
