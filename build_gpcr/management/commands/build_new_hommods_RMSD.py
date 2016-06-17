# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 14:14:58 2016

@author: gaspar
"""

from django.core.management.base import BaseCommand

import structure.list_new_unique_xtals as nx
from build_gpcr.management.commands.build_homology_models import *
from structure.management.commands.calculate_RMSD import *
from structure.models import StructureModel
import os
import Bio.PDB as PDB
import requests
import zipfile
import numpy as np
import re


class Command(BaseCommand):
    def handle(self, *args, **options):
        v = Validation()
        v.run_RMSD_list(['./structure/PDB/pdb4zjc.ent','./structure/homology_models/ox1r_human_Inactive/ox1r_human_Inactive_model_updated.pdb','./structure/homology_models/SwissMod/ox1r_human_2016-05-04/model/01/ox1r_model.pdb'],[99, 102, 103, 112, 122, 123, 126, 130, 179, 183, 204, 219, 314, 318, 344, 348])
#        v.run_RMSD_list(['./structure/PDB/pdb5cxv.ent','./structure/homology_models/acm1_human_Inactive/acm1_human_Inactive_model.pdb','./structure/homology_models/SwissMod/acm1_human_2016-05-27/model/01/acm1_model.pdb','./structure/homology_models/SwissMod/acm1_human_2016-05-27/model/acm1_model.pdb'],[105, 106, 109, 110, 157, 189, 192, 193, 196, 197, 378, 381, 382, 404, 407, 408])
#        files = []        
#        gpcrm_dir = os.listdir('./structure/homology_models/GPCRM/ACM1_HUMAN_inact')
#        for i in gpcrm_dir:
#            if i.startswith('gpcrm'):
#                files.append(os.path.join('./structure/homology_models/GPCRM/ACM1_HUMAN_inact',i))        
#        v.run_RMSD_list(['./structure/PDB/pdb5dsg.ent','./structure/homology_models/acm4_human_Inactive/acm4_human_Inactive_model.pdb','./structure/homology_models/SwissMod/acm4_human_2016-05-27/model/01/acm4_model.pdb','./structure/homology_models/SwissMod/acm4_human_2016-05-27/model/acm4_model.pdb'],[112, 113, 116, 117, 164, 190, 196, 199, 200, 203, 204, 413, 416, 417, 439, 442, 443])
        c=0
        ar = np.array([0,0,0,0])
        for k,l in v.number_of_residues_superposed.items():
            print(k)
            print(l)
        for i,j in v.rmsds.items():
            c+=1
            if c>4:
#                print([j['overall_all'],j['overall_backbone'],j['TM_all'],j['TM_backbone']])
                ar = np.vstack((ar,[j['overall_all'],j['overall_backbone'],j['TM_all'],j['TM_backbone']]))
            else:
                print(i)
                pprint.pprint(j)
        ar = ar[1:]
        print(ar)
        print(np.mean(ar,axis=0))
#        n = nx.QueryPDB()
#        n.list_xtals()
#        for r in n.new_uniques:
#            try:
#                hm = list(StructureModel.objects.filter(protein__entry_name=r[1].entry_name, state__name="Inactive"))[-1]
#                with open('./structure/homology_models/model.pdb', 'w+') as f:
#                    f.write(hm.pdb)
#                gpcrdb_path = './structure/homology_models/model.pdb'
#                print('Fetched from db')
#            except:
#                if os.path.isfile('./structure/homology_models/'+r[1].entry_name+'_Inactive/'+r[1].entry_name+'_Inactive_model.pdb'):               
#                    gpcrdb_path = './structure/homology_models/{}_Inactive/{}_Inactive_model.pdb'.format(r[1].entry_name,r[1].entry_name)
#                else:
#                    try:
#                        hm = HomologyModeling(r[1].entry_name,'Inactive',['Inactive'])
#                        alignment = hm.run_alignment()
#                        hm.build_homology_model(alignment)
#                        hm.format_final_model()
#                        gpcrdb_path = './structure/homology_models/{}_Inactive/{}_Inactive_model.pdb'.format(r[1].entry_name,r[1].entry_name)
#                    except:
#                        print('ERROR with model {}'.format(r))
#                        gpcrdb_path = None
#                        continue
#            with open('./static/homology_models/RMSD_table.csv', 'r') as csvfile:
#                lines = csvfile.readlines()
#                present = False
#                for line in lines:
#                    if r[0] in line and 'v{}'.format(hm.version) in line:
#                        present = True
#            if present==True:
#                continue
#            if not os.path.exists('./structure/PDB/pdb{}.ent'.format(r[0].lower())):
#                l = PDB.PDBList()
#                l.retrieve_pdb_file(r[0],pdir='./structure/PDB/')
#
#            gpcrm_link = 'http://gpcrm.biomodellab.eu/models/data/Models_of_{}_receptors/{}_{}.zip'.format('inactive',r[1].entry_name.upper(),'inact') 
#            gpcrm_path = './structure/homology_models/GPCRM/{}_inact.zip'.format(r[1].entry_name.upper())
#            try:
#                if not os.path.isfile(gpcrm_path):
#                    req = requests.get(gpcrm_link)
#                    if req.status_code!=200:
#                        raise Exception()
#                    with open(gpcrm_path, "wb") as code:
#                        code.write(req.content)
#                z = zipfile.ZipFile(gpcrm_path)
#                new_path = os.path.join(gpcrm_path.split('.zip')[0])
#                z.extractall(new_path)
#                files = os.listdir(new_path)
#                file_count = 0
#                ar = np.array([0.0,0.0,0.0,0.0])
#                for f in files:
#                    if f.startswith('gpcrm_model'):
#                        file_count+=1
#                        v1 = Validation()
#                        rmsds = v1.run_RMSD('./structure/PDB/pdb{}.ent'.format(r[0].lower()),new_path+'/'+f)
#                        ar+=np.array(rmsds)
#                best_result = ar/file_count
#                print(r)
#                print('GPCRM (averaged): ',best_result)
#            except:
#                pass
#            if gpcrdb_path!=None:
#                v2 = Validation()
#                rmsds = v2.run_RMSD('./structure/PDB/pdb{}.ent'.format(r[0].lower()),
#                                    gpcrdb_path)
#                if gpcrdb_path.split('/')[-1]=='model.pdb':
#                    os.remove(gpcrdb_path)
#                print('GPCRdb: ',rmsds)
#            date = re.search('\d{4}/\d{2}/\d{2}',hm.pdb).group(0)
#            rmsds = [r[1].entry_name, '{}_GPCRdb'.format(r[0])] + [str('%.1f' % round(i, 1)) for i in rmsds] + ['v{}'.format(hm.version), date]
#            best_result = [r[1].entry_name, '{}_GPCRM'.format(r[0])] + [str('%.1f' % round(i, 1)) for i in best_result] + ['v{}'.format(hm.version), 'None']
#            swiss = [r[1].entry_name, '{}_SwissModel'.format(r[0]),'None','None','None','None'] + ['v{}'.format(hm.version), 'None']
#            
#            with open('./static/homology_models/RMSD_table.csv', 'a') as csvfile:
#                csvfile.write(', '.join(rmsds)+'\n')
#                csvfile.write(', '.join(best_result)+'\n')
#                csvfile.write(', '.join(swiss)+'\n')
#                print('CSV file done')

#        v = Validation()
#        rmsds = v.run_RMSD('./structure/homology_models/ox1r_human_Inactive/ox1r_human_Inactive_model.pdb',
#                           './structure/PDB/pdb4zj8.ent')
#        print('GPCRdb: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/GoMoDo/ox1r_human_Inactive/ox1r_human.B99990001_2z73.pdb',
#                           './structure/PDB/pdb4zj8.ent')
#        print('GoMoDo 2z73: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/GoMoDo/ox1r_human_Inactive/ox1r_human.B99990001_4iar.pdb',
#                           './structure/PDB/pdb4zj8.ent')
#        print('GoMoDo 4iar: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/SwissMod/ox1r_human_2016-05-04/model/01/ox1r_model.pdb',
#                           './structure/PDB/pdb4zj8.ent')
#        print('SwissModel 4zj8: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/SwissMod/ox1r_human_2016-05-04/model/01/ox1r_model.pdb',
#                           './structure/PDB/pdb4zjc.ent')
#        print('SwissModel 4zjc: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/SwissMod/acm1_human_2016-05-27/model/01/acm1_model.pdb',
#                           './structure/PDB/pdb5cxv.ent')
#        print('SwissModel 5CXV: ',rmsds)
#        rmsds = v.run_RMSD('./structure/homology_models/SwissMod/acm4_human_2016-05-27/model/01/acm4_model.pdb',
#                           './structure/PDB/pdb5dsg.ent')
#        print('SwissModel 5DSG: ',rmsds)
        


        