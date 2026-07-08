from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from django.conf import settings

from django.db.models import Count, Max, F, Subquery, OuterRef

from  sqlite3 import IntegrityError

from ligand.models import Ligand

import os
import django.apps
import logging
import hashlib
import base64
import sqlite3
import csv
import copy

TABLE_NAME = 'ligand' #cannot be "U0"
db_file_path_default = os.path.join(settings.DATA_DIR,"structure_data","ligand_sequence_hash.sqlite3")
max_batch_size = 1000


def _insert_into_table_bulk_sqlite(connection,cursor,table_name,update_fields,qvalues):
    sql_query_insert_into = "INSERT INTO %s (%s) values (%s)" % (table_name,','.join([f for f in update_fields]),','.join(['?']*len(update_fields)))
    cursor.executemany(sql_query_insert_into,[[qvalue[f] for f in update_fields] for qvalue in qvalues])
    connection.commit()

def _update_table_bulk_sqlite(connection,cursor,table_name,query_fields,update_fields,qvalues):
    sql_query_insert_into = "UPDATE %s SET %s WHERE %s" % (table_name,', '.join([f+' = ?' for f in update_fields]),', '.join([f+' = ?' for f in query_fields]))
    cursor.executemany(sql_query_insert_into,[[qvalue[f] for f in update_fields+query_fields] for qvalue in qvalues])
    connection.commit()

def _initialize_db_table(con):
        cur = con.cursor()
        cur.execute(
            "DROP TABLE IF EXISTS %s"  % (TABLE_NAME)
        )
        cur.close()
        cur = con.cursor()
        cur.execute(
            "CREATE TABLE IF NOT EXISTS %s(" % (TABLE_NAME) + \
                "id INTEGER PRIMARY KEY, "
                "gpcrdb_pk INTEGER UNIQUE, "
                "sequence TEXT NOT NULL, "
                "sequence_hash TEXT NOT NULL, "
                "sequence_hash_col TEXT DEFAULT \"0\", "
                "sequence_dup INTEGER NOT NULL DEFAULT 0, "
                "sequence_hash_and_col_main BOOLEAN NOT NULL DEFAULT 1 CHECK (sequence_hash_and_col_main IN (0, 1)),"
                "UNIQUE (sequence_hash, sequence_hash_col, sequence_dup)"
            ")"
        )# "UNIQUE (sequence_hash, sequence_hash_col, sequence_dup)"
        cur.close()

        cur = con.cursor()
        cur.execute("CREATE INDEX IF NOT EXISTS sequence_hash_index ON %s (sequence_hash)" % (TABLE_NAME))
        cur.execute("CREATE INDEX IF NOT EXISTS sequence_hash_sequence_hash_col_index ON %s (sequence_hash,sequence_hash_col)" % (TABLE_NAME))
        cur.close()

def _create_artificial_collision(q_results, current_batch_start):
    col_i = 0
    q_result_n_1 = copy.deepcopy(q_results[-1])
    
    if current_batch_start == 0:
        collision_test_hash = q_result_n_1['sequence_hash'] = q_results[-2]['sequence_hash']
        col_sequence = q_result_n_1['sequence']
    else:
        # make sure sequence is different
        if col_sequence == q_result_n_1['sequence']:
            for q_result in q_results:
                if  q_result['sequence'] != col_sequence:
                    col_sequence = q_result_n_1['sequence'] = q_result['sequence']

        q_result_n_1['sequence_hash'] = collision_test_hash
    q_result_n_1['gpcrdb_pk'] = 2147483647 - col_i
    col_i += 1
    q_results.append(q_result_n_1)
    for j in range(0,3):
        q_result_n_1_dup = copy.deepcopy(q_result_n_1)
        q_result_n_1_dup['gpcrdb_pk'] = 2147483647 - col_i
        col_i += 1
        q_results.append(q_result_n_1_dup)
    del q_result_n_1_dup
    del q_result_n_1
    print('Test collision on hash: '+collision_test_hash)  

def _result_set_to_hash_keyed_dict(query_result_set):
    """Create a dictionary that uses sequence hashes as keys and a list of Ligand objects with the same hash as values"""
    sequence_hashes_dict = {}
    for q_result in query_result_set:
        my_hash = q_result['sequence_hash']
        if my_hash in sequence_hashes_dict:            
            sequence_hashes_dict[my_hash].append(q_result)
        else:
            sequence_hashes_dict[my_hash] = [q_result]
    return sequence_hashes_dict

class Command(BaseCommand):
    help = 'Build ligand sequence hash. '

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--output', default=db_file_path_default, action='store', help='Output path of the sqlite3 file.')
        parser.add_argument('--verbose', default=False, action='store_true', help='Print progress in stdout.')
        parser.add_argument('--collision-test', default=False, action='store_true', help='Only for code testing.')
        parser.add_argument('--debug-csv', default=False, action='store_true', help='Creates a debug CSV file.')

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):    
        error = None
        db_file_path = options['output']
        if options['verbose']: print('Building ligand sequence hashes...')
        
        con = sqlite3.connect(db_file_path)

        _initialize_db_table(con)

        query_fields = ['id','sequence']
        update_fields = ['gpcrdb_pk','sequence','sequence_hash','sequence_hash_col','sequence_dup','sequence_hash_and_col_main']
        q = Ligand.objects.exclude(sequence=None).order_by('id').values(*query_fields)
        current_batch_start = 0        
        collision_test_hash = None
        col_sequence = None
        while True: # Run in batches of batch_size to be memory efficient
            q_results = list(q[current_batch_start:current_batch_start + max_batch_size])
            if not q_results:
                break
            if options['verbose']: print('Parsing from '+str(current_batch_start+1)+' to '+str(current_batch_start+len(q_results)))
            for q_result in q_results:    
                my_hash = base64.b32encode(hashlib.md5(q_result['sequence'].encode()).digest()).decode().strip('=')
                q_result['gpcrdb_pk'] = q_result['id']
                q_result['id'] = None
                q_result['sequence_hash'] = my_hash
                q_result['sequence_hash_col'] = '0'
                q_result['sequence_dup'] = 0
                q_result['sequence_hash_and_col_main'] = 1

            # collision test for development
            if options['collision_test']:
                _create_artificial_collision(q_results, current_batch_start)
            
            # Try to save the hashes in SQLlite DB 
            cur = con.cursor()
            try:
                _insert_into_table_bulk_sqlite(con,cur,TABLE_NAME,update_fields,q_results)
                cur.close()

            except IntegrityError as e:
                # This runs on duplicates or collisions
                con.rollback()
                cur.close()                
                
                sequence_hashes_dict = _result_set_to_hash_keyed_dict(q_results)
                
                #extract the sequence hashes that have duplicates/collisions
                dup_sequence_hashes_set = set([key for key, value in sequence_hashes_dict.items() if len(value) > 1])

                list_of_unique_hashes = list(sequence_hashes_dict.keys())
                q_num_col = con.cursor()

                # This SQL query returns a table with sequence hashes, hash collision ID and the number of duplicates 
                sql_query = 'SELECT "table_name_0"."sequence_hash", COUNT("table_name_0"."id") AS "num_col"'+ \
                            ' FROM "%s" AS "table_name_0"' % (TABLE_NAME)  + \
                            ' WHERE "table_name_0"."sequence_dup" = (' + \
                                'SELECT MAX(U0."sequence_dup") AS "max_sequence_dup" ' + \
                                'FROM "%s" U0 ' % (TABLE_NAME) + \
                                'WHERE (' + \
                                'U0."sequence_hash" IN ('+','.join(['?']*len(list_of_unique_hashes))+') ' + \
                                'AND U0."sequence_hash" = ("table_name_0"."sequence_hash") ' + \
                                    'AND U0."sequence_hash_col" = ("table_name_0"."sequence_hash_col")' + \
                                ') '+ \
                                'GROUP BY U0."sequence_hash", U0."sequence_hash_col"' + \
                            ') GROUP BY "table_name_0"."sequence_hash"'
                q_num_col.execute(sql_query,list_of_unique_hashes)
                
                # Make a list of batches of sequence hashes than will return a query with a number of records < batch_size
                # The following code does not preserve the order of the sequence hashes
                col_batch_list = []
                overrun_buffer_hashes_list = []
                current_batch_size = 0
                batch = []
                for row in q_num_col:
                    num_col = row[1]
                    sequence_hash = row[0]
                    if num_col > max_batch_size:
                        
                        # Remove this warning if the code already takes care of this
                        msg = "build_ligand_sequence_hash: Number of " + \
                                            "hash collisions larger than %d for %s." % (max_batch_size, sequence_hash) + \
                                            "This might cause RAM memory overrun during building with the current " + \
                                            "implementation."
                        self.logger.warning(msg)
                        if options['verbose']: print('WARNING: '+msg)
                        overrun_buffer_hashes_list.append(sequence_hash)
                        continue
                    
                    current_batch_size += num_col
                    if current_batch_size > max_batch_size:
                        col_batch_list.append(batch)
                        batch = []
                        current_batch_size = num_col
                    batch.append(sequence_hash)
                if len(batch) > 0:
                    col_batch_list.append(batch)

                # Today, we still don't take care of buffer overruning hashes
                col_batch_list += overrun_buffer_hashes_list
                del overrun_buffer_hashes_list
                q_num_col.close()
                del q_num_col
                
                fixed_dup_or_col_hashes = set()
                records_to_update = []
                for batch in col_batch_list:
                    
                    sql_query = 'SELECT %s ' % (','.join(['"table_name_0"."'+f+'"' for f in ['id']+update_fields]))+ \
                            'FROM "%s" AS "table_name_0"' % (TABLE_NAME)  + \
                            'WHERE "table_name_0"."sequence_dup" = (' + \
                                'SELECT MAX(U0."sequence_dup") AS "max_sequence_dup" ' + \
                                'FROM "%s" U0 ' % (TABLE_NAME) + \
                                'WHERE (' + \
                                    'U0."sequence_hash" IN ('+','.join(['?']*len(batch))+') ' + \
                                    'AND U0."sequence_hash" = ("table_name_0"."sequence_hash") ' + \
                                    'AND U0."sequence_hash_col" = ("table_name_0"."sequence_hash_col")' + \
                                ') '+ \
                                'GROUP BY U0."sequence_hash", U0."sequence_hash_col"' + \
                            ')'

                    q_col = con.execute(sql_query,batch)

                    q_col_or_dup_hash_seq_dict = {}
                    for raw_r in q_col:
                        fields = ['id']+update_fields
                        r = { fields[j] : r_val for j, r_val in enumerate(raw_r)}
                        my_hash = r['sequence_hash']
                        if my_hash not in q_col_or_dup_hash_seq_dict:
                            q_col_or_dup_hash_seq_dict[my_hash] = {}
                        q_col_or_dup_hash_seq_dict[my_hash][r['sequence']] = r
                    del q_col

                    for my_hash, seq_dict in q_col_or_dup_hash_seq_dict.items():
                        # hash_cols is an hexadecimal number written from right to left
                        hash_cols_decimal_max = max([int(r['sequence_hash_col'][::-1],16) for r in seq_dict.values()]\
                                                      +[0])
                        
                            
                        db_col_dup_count_dict = {}  # duplicate count for sequences already in SQLlite DB
                        new_col_dup_count_dict = {} # duplicate count for new collisions not in SQLlite DB
                        new_col_count = 0           # hash collision count for new collisions not in SQLlite DB
                        new_col_sequences_dict = {} # sequences that have new hash collisions not in SQLlite DB
                        sequences = list(seq_dict.keys())
                        sequences_num = len(sequences)
                        for q_result in sequence_hashes_dict[my_hash]:
                            # check if it is a duplicate
                            if  q_result['sequence'] in seq_dict:
                                # It is a duplicate of a hash already stored in the SQLlite DB
                                r = seq_dict[q_result['sequence']]
                                hash_col = r['sequence_hash_col']
                                q_result['sequence_hash_col'] = hash_col
                                if hash_col not in db_col_dup_count_dict:
                                    db_col_dup_count_dict[hash_col] = r['sequence_dup']
                                db_col_dup_count_dict[hash_col] += 1
                                q_result['sequence_dup'] = db_col_dup_count_dict[hash_col]
                                q_result['sequence_hash_and_col_main'] = False
                                


                            elif q_result['sequence'] in new_col_sequences_dict:
                                # It is a duplicate of a hash of a new hash collisions not in SQLlite DB
                                hash_col = new_col_sequences_dict[q_result['sequence']]
                                q_result['sequence_hash_col'] = hash_col
                                if new_col_count not in new_col_dup_count_dict:
                                    new_col_dup_count_dict[new_col_count] = 0
                                new_col_dup_count_dict[new_col_count] += 1
                                q_result['sequence_dup'] = new_col_dup_count_dict[new_col_count]
                                q_result['sequence_hash_and_col_main'] = False
                            else:
                                # It is a collision
                                new_col_count += 1
                                hash_col = hex(hash_cols_decimal_max + new_col_count)[2:][::-1]
                                q_result['sequence_hash_col'] = hash_col
                                q_result['sequence_hash_and_col_main'] = True
                                q_result['sequence_dup'] = 0
                                new_col_sequences_dict[q_result['sequence']] = hash_col
                        if hash_cols_decimal_max == 0:
                            # Update, in SQL lite DB, hash colision IDs values that were empty
                            if sequences_num > 1:
                                msg = "build_ligand_sequence_hash: sequence hash %s " % (sequence_hash) + \
                                      "had collisions and sequence_hash_col field is empty."
                                self.logger.warning(msg)
                                if options['verbose']: print('WARNING: '+msg)

                            if new_col_count > 0:
                                q_result = seq_dict[sequences[0]]
                                q_result['sequence_hash_col'] = '0'
                                records_to_update.append(q_result)
                            
                        fixed_dup_or_col_hashes.add(my_hash)
                        
                cur = con.cursor()
                _update_table_bulk_sqlite(con,cur,TABLE_NAME,['id'],update_fields,records_to_update)
                cur.close()
                dup_or_col_hashes_not_in_db = dup_sequence_hashes_set - fixed_dup_or_col_hashes  
                
                # Assign duplicated hashes/collisions IDs for hashes no it the SQL lite DB
                for my_hash in list(dup_or_col_hashes_not_in_db):
                    dup_col_objs = sequence_hashes_dict[my_hash]
                    new_col_dup_count_dict = {0:0}
                    new_col_count = 0
                    new_col_sequences_dict = {}
                    for q_result in dup_col_objs:
                        if q_result['sequence'] in new_col_sequences_dict:
                            # It is a duplicate
                            hash_col = new_col_sequences_dict[q_result['sequence']]
                            q_result['sequence_hash_col'] = hash_col
                            new_col_dup_count_dict[new_col_count] += 1
                            q_result['sequence_dup'] = new_col_dup_count_dict[new_col_count]
                            q_result['sequence_hash_and_col_main'] = False
                        else:
                            # It is a collision
                            hash_col = hex(new_col_count)[2:][::-1]
                            q_result['sequence_hash_col'] = hash_col
                            q_result['sequence_hash_and_col_main'] = True
                            q_result['sequence_dup'] = 0
                            new_col_sequences_dict[q_result['sequence']] = hash_col
                            new_col_count += 1
                            new_col_dup_count_dict[new_col_count] = 0
                    if len(new_col_sequences_dict.keys()) < 2:
                        # Update hash colision IDs values that were empty
                        for q_result in dup_col_objs:
                            q_result['sequence_hash_col'] = '0'
                            
                # Convert into a list the dict() with a list of Ligand objects with duplicated hashes/collisions
                q_results = []
                for v in sequence_hashes_dict.values():
                    q_results += v
                       
                del sequence_hashes_dict
                if options['debug_csv']:
                    with open(db_file_path+'_'+str(current_batch_start)+'.csv', 'w', newline='') as csvfile:
                        fieldnames = ['id']+update_fields
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for r in q_results:
                            writer.writerow(r)
                cur = con.cursor()
                _insert_into_table_bulk_sqlite(con,cur,TABLE_NAME,update_fields,q_results)
                cur.close()



            current_batch_start += max_batch_size
            
        con.close()
        self.logger.info('Ligand sequence hashes built.')
            



 