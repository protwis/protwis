
import xlrd

import json
import urllib.request
#env/bin/python3 -m pip install xlrd

from mutation.models import *
from datetime import datetime

def loaddatafromexcel(excelpath):
	workbook = xlrd.open_workbook(excelpath)
	worksheets = workbook.sheet_names()
	temp = []
	for worksheet_name in worksheets:
		worksheet = workbook.sheet_by_name(worksheet_name)
		#print(worksheet_name)

		if worksheet.cell_value(0, 0) == "REFERENCE \nDOI (or PMID)":
			#print("MATCH")
			pass
		else:
			#print("SKIPPING")
			continue

		num_rows = worksheet.nrows - 1
		num_cells = worksheet.ncols - 1
		curr_row = 0 #skip first, otherwise -1
		while curr_row < num_rows:
			curr_row += 1
			row = worksheet.row(curr_row)
			#print('Row:', curr_row)
			curr_cell = -1
			temprow = []
			if worksheet.cell_value(curr_row, 0) == '': #if empty
				continue
			while curr_cell < num_cells:
				curr_cell += 1
				# Cell Types: 0=Empty, 1=Text, 2=Number, 3=Date, 4=Boolean, 5=Error, 6=Blank
				cell_type = worksheet.cell_type(curr_row, curr_cell)
				cell_value = worksheet.cell_value(curr_row, curr_cell)
				#print('	', cell_type, ':', cell_value)
				temprow.append(cell_value)
			temp.append(temprow)
			#if curr_row>10: break
		return temp

def analyse_rows(rows):
	temp = []
	for r in rows:
		d = {}
		d['reference'] = r[0]
		d['protein'] = r[1].replace("__","_").lower()
		d['mutation_pos'] = r[2]
		d['mutation_from'] = r[3]
		d['mutation_to'] = r[4]
		d['ligand_name'] = r[5]
		d['ligand_type'] = r[6]
		d['ligand_id'] = r[7]
		d['ligand_class'] = r[8]
		d['exp_type'] = r[9]
		d['exp_func'] = r[10]
		d['exp_wt_value'] = int(r[11]) if r[11] else 0
		d['exp_wt_unit'] = r[12]
		d['exp_mu_effect_type'] = r[13]
		d['exp_mu_effect_sign'] = r[14]
		d['exp_mu_value_raw'] = int(r[15]) if r[15] else 0
		d['exp_mu_effect_qual'] = r[16]
		d['exp_mu_effect_ligand_prop'] = r[17]
		d['exp_mu_ligand_ref'] = r[18]
		d['opt_type'] = r[21]
		d['opt_wt'] = int(r[22]) if r[22] else 0
		d['opt_mu'] = int(r[23]) if r[23] else 0
		d['opt_sign'] = r[24]
		d['opt_percentage'] = int(r[25]) if r[25] else 0
		d['opt_qual'] = r[26]
		d['opt_agonist'] = r[27]



		if isinstance(d['ligand_id'], float): d['ligand_id'] = int(d['ligand_id'])
		if isinstance(d['mutation_pos'], float): d['mutation_pos'] = int(d['mutation_pos'])


		temp.append(d)
	return temp


def insert_raw(r):
	obj, created = MutationRaw.objects.get_or_create(
	reference=r['reference'], 
	protein=r['protein'], 
	mutation_pos=r['mutation_pos'], 
	mutation_from=r['mutation_from'], 
	mutation_to=r['mutation_to'], 
	ligand_name=r['ligand_name'], 
	ligand_idtype=r['ligand_type'], 
	ligand_id=r['ligand_id'], 
	ligand_class=r['ligand_class'], 
	exp_type=r['exp_type'], 
	exp_func=r['exp_func'], 
	exp_wt_value=r['exp_wt_value'], #
	exp_wt_unit=r['exp_wt_unit'], 
	exp_mu_effect_type=r['exp_mu_effect_type'], 
	exp_mu_effect_sign=r['exp_mu_effect_sign'], 
	exp_mu_effect_value=r['exp_mu_value_raw'], #
	exp_mu_effect_qual=r['exp_mu_effect_qual'], 
	exp_mu_effect_ligand_prop=r['exp_mu_effect_ligand_prop'], 
	exp_mu_ligand_ref=r['exp_mu_ligand_ref'], 
	opt_type=r['opt_type'], 
	opt_wt=r['opt_wt'], #
	opt_mu=r['opt_mu'], #
	opt_sign=r['opt_sign'], 
	opt_percentage=r['opt_percentage'], #
	opt_qual=r['opt_qual'], 
	opt_agonist=r['opt_agonist'], 
	added_by='munk', 
	added_date=datetime.now()
	)

	raw_id = obj

	return raw_id


def check_reference(r):
	
	ref=MutationRefs.objects.filter(reference=r)


	if ref.exists():
		ref=MutationRefs.objects.get(reference=r)
		#return ref.id
		return ref
	else:
		url = "http://search.crossref.org/dois?q="+r
		
		response = urllib.request.urlopen(url)

		j = json.loads(response.read().decode('utf-8'))

		if len(j)>0:
			j = j[0]
			e = MutationRefs(ref_type='DOI', link=j['doi'], title=j['title'], citation= j['fullCitation'], year= j['year'], reference=r)
			e.save()

			#return e.id
			return e
		else:
			print('No reference?')
			print(r)
			return None

def get_ligand(r):
	
	check=MutationLigand.objects.filter(idtype=r['ligand_type'], idid=r['ligand_id'])


	if check.exists():
		check=MutationLigand.objects.get(idtype=r['ligand_type'], idid=r['ligand_id'])
		return check
	else:

		if r['ligand_type'] == 'PubChem CID':
			#https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/57328469/property/IUPACName,CanonicalSMILES,IsomericSMILES/JSON
			url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(int(r['ligand_id']))+"/property/IUPACName,CanonicalSMILES,IsomericSMILES/JSON";
			response = urllib.request.urlopen(url)
			j = json.loads(response.read().decode('utf-8'))

			if j['PropertyTable']['Properties'][0]['CID']:
				smiles = j['PropertyTable']['Properties'][0]['CanonicalSMILES'];
			else:
				smiles = r['ligand_id']

			e = MutationLigand(idtype=r['ligand_type'], idid=r['ligand_id'], name=r['ligand_name'], longseq= smiles)
			e.save()
			return e
		elif r['ligand_type'] == 'SMILES':
			e = MutationLigand(idtype=r['ligand_type'], idid=r['ligand_id'], name=r['ligand_name'], longseq=r['ligand_id'])
			e.save()
			return e
		else:
			return None


		#$q = "INSERT INTO ligands (idtype,name,idid,longseq) VALUES ('$ligand_type','$ligand_name','$ligand_id','$smiles')";

		

