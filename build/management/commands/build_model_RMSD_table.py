# -*- coding: utf-8 -*-
"""
Created on Sat May 28 17:36:45 2016

@author: Gaspar Pandy
"""

from build.management.commands.base_build import Command as BaseBuild

from structure.models import StructureModel, StructureModelRMSD
import datetime

class Command(BaseBuild):
    def handle(self, *args, **options):
        with open('./static/homology_models/RMSD_table.csv', 'r') as csvfile:
            lines = csvfile.readlines()
            for line in lines[1:]:
                data = line.split(',')
                service = data[1].split('_')[-1]
                if service=='GPCRdb':
                    hm = StructureModel.objects.get(protein__entry_name=data[0], state__name='Inactive', version=float(data[6][1:]))
                try:
                    StructureModelRMSD.get(homology_model=hm, service=service)
                except:
                    d = data[7].split('/')
                    d = datetime.date(int(d[0]),int(d[1]),int(d[2]))
                    rmsd, created = StructureModelRMSD.objects.update_or_create(service=service, homology_model=hm,
                                                                                pdb=data[1].split('_')[0],
                                                                                overall_all=float(data[2]), 
                                                                                overall_backbone=float(data[3]),
                                                                                TM_all=float(data[4]),
                                                                                TM_backbone=float(data[5]),
                                                                                version=float(data[6][1:]),
                                                                                date = d)
                    print('{} not in table'.format(data))