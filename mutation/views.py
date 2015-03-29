from django.shortcuts import render

from django.http import HttpResponse
from mutation.testing import *

from mutation.models import *

from residue.models import Residue
from protein.models import Protein


from datetime import datetime
#env/bin/python3 -m pip install xlrd

import re
import math


# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")

def importmutation(request):

	rows = loaddatafromexcel('/vagrant/shared/protwis/mutation/import.xlsx')

	rows = analyse_rows(rows)

	whattoreturn = []
	c = 0
	skipped = 0
	inserted = 0
	for r in rows:
		raw_id = insert_raw(r)
		ref_id = check_reference(r['reference'])
		lig_id = get_ligand(r)

		protein_id = 0
		residue_id = 0

		check=Protein.objects.filter(entry_name=r['protein'])
		if check.exists():
			check=Protein.objects.get(entry_name=r['protein'])
			protein_id = check
		else:
			whattoreturn.append(['Skipped due to no protein',r['protein']])
			skipped += 1
			continue

		check=Residue.objects.filter(protein=protein_id,sequence_number=r['mutation_pos'])
		if check.exists():
			check=Residue.objects.get(protein=protein_id,sequence_number=r['mutation_pos'])
			residue_id = check
		else:
			whattoreturn.append(['Skipped due to no residue',r['protein'],r['mutation_pos']])
			skipped += 1
			continue
			# residue_id = Residue()
			# residue_id.protein = protein_id
			# residue_id.sequence_number = r['mutation_pos']
			# residue_id.amino_acid = r['mutation_from']  
			# residue_id.save()

		obj, created = MutationLigandClass.objects.get_or_create(classname=r['ligand_class'])
		ligclass_id = obj


		obj, created = MutationLigandRef.objects.get_or_create(reference=r['exp_mu_ligand_ref'])
		ligref_id = obj

		obj, created = MutationType.objects.get_or_create(type=r['exp_type'])
		exp_type_id = obj

		obj, created = MutationFunc.objects.get_or_create(func=r['exp_func'])
		exp_func_id = obj

		obj, created = MutationMeasure.objects.get_or_create(measure=r['exp_mu_effect_type'])
		exp_measure_id = obj

		obj, created = MutationQual.objects.get_or_create(qual=r['exp_mu_effect_qual'], prop=r['exp_mu_effect_ligand_prop'])
		exp_qual_id = obj

		obj, created = MutationLigandRef.objects.get_or_create(reference=r['exp_mu_ligand_ref'])
		effect_ligand_reference_id = obj

		obj, created = 	MutationOptional.objects.get_or_create(type=r['opt_type'], wt=r['opt_wt'], mu=r['opt_mu'], sign=r['opt_sign'], percentage=r['opt_percentage'], qual=r['opt_qual'], agonist=r['opt_agonist'])
		exp_opt_id = obj


		
		logtypes = ['pEC50','pIC50','pK']
		
		
		foldchange = 0
		typefold = ''
		if r['exp_mu_effect_type']=='Activity/affinity' and r['exp_wt_unit']!=0:
					
			if re.match("(" + ")|(".join(logtypes) + ")", r['exp_type']):  #-log values!
				foldchange = round(math.pow(10,-r['exp_mu_value_raw'])/pow(10,-r['exp_wt_value']),3);
				typefold = r['exp_type']+"_log"
			else:
				foldchange = round(r['exp_mu_value_raw']/r['exp_wt_value'],3);
				typefold = r['exp_type']+"_not_log"
			
			
			if foldchange<1 and foldchange!=0:
				foldchange = -round((1/foldchange),3)
			elif r['exp_mu_effect_type'] =='Fold effect (mut/wt)':
				foldchange = round(r['exp_mu_value_raw'],2);
				if foldchange<1: foldchange = -round((1/foldchange),3);
		

		obj, created = Mutation.objects.get_or_create(
		refs=ref_id, 
		protein=protein_id, 
		residue=residue_id, #MISSING 
		ligand=lig_id, 
		ligand_class=ligclass_id, 
		ligand_ref = ligref_id,
		raw = raw_id,
		optional = exp_opt_id,
		exp_type=exp_type_id, 
		exp_func=exp_func_id, 
		exp_measure = exp_measure_id,
		exp_qual = exp_qual_id,


		mutation_to=r['mutation_to'], 
		wt_value=r['exp_wt_value'], #
		wt_unit=r['exp_wt_unit'], 

		mu_value = r['exp_mu_value_raw'],
		mu_sign = r['exp_mu_effect_sign'], 
		foldchange = foldchange

		#added_by='munk', 
		#added_date=datetime.now()
		)

		mut_id = obj.id


		#whattoreturn.append([protein_id,residue_id,raw_id,ref_id,lig_id,ligclass_id,exp_type_id,exp_func_id,exp_measure_id,exp_qual_id,typefold,foldchange,mut_id])
		whattoreturn.append(['Inserted',protein_id.entry_name,residue_id.sequence_number,foldchange])
		inserted += 1
		c += 1
		#if c>10: break


	context = {'rows': whattoreturn, 'skipped' : skipped, 'inserted' : inserted}

	#return HttpResponse(ref_id)
	return render(request,'mutation/index.html',context)