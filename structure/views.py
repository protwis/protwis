from django.shortcuts import render
from django.conf import settings
from django.views.generic import TemplateView, View
from django.http import HttpResponse, HttpResponseRedirect
from django.db.models import Count, Q, Prefetch, TextField
from django.db.models.functions import Concat
from django import forms

from django.shortcuts import redirect

from common.phylogenetic_tree import PhylogeneticTreeGenerator
from protein.models import ProteinSegment
from structure.models import Structure, StructureModel, StructureComplexModel, StructureExtraProteins, StructureVectors, StructureModelRMSD
from structure.functions import CASelector, SelectionParser, GenericNumbersSelector, SubstructureSelector, ModelRotamer
from structure.assign_generic_numbers_gpcr import GenericNumbering, GenericNumberingFromDB
from structure.structural_superposition import ProteinSuperpose,FragmentSuperpose
from structure.forms import *
from signprot.models import SignprotComplex, SignprotStructure, SignprotStructureExtraProteins
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from protein.models import Protein, ProteinFamily
from construct.models import Construct
from construct.functions import convert_ordered_to_disordered_annotation,add_construct
from common.views import AbsSegmentSelection,AbsReferenceSelection
from common.selection import Selection, SelectionItem
from common.extensions import MultiFileField
from common.models import ReleaseNotes
from common.alignment import Alignment, GProteinAlignment
from residue.models import Residue

import io
import numpy as np
from Bio.PDB.vectors import Vector, rotmat

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

import inspect
import os
import time
import zipfile
import json

from copy import deepcopy
from io import StringIO, BytesIO
from collections import OrderedDict
from Bio.PDB import PDBIO, PDBParser

from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
import smtplib

class_dict = {'001':'A','002':'B1','003':'B2','004':'C','005':'D1','006':'F','007':'T','008':'O'}

class StructureBrowser(TemplateView):
	"""
	Fetching Structure data for browser
	"""
	template_name = "structure_browser.html"

	def get_context_data (self, **kwargs):

		context = super(StructureBrowser, self).get_context_data(**kwargs)
		try:
			structures = Structure.objects.all().select_related(
				"state",
				"pdb_code__web_resource",
				"protein_conformation__protein__species",
				"protein_conformation__protein__source",
				"protein_conformation__protein__family__parent__parent__parent",
				"publication__web_link__web_resource").prefetch_related(
				"stabilizing_agents", "construct__crystallization__crystal_method",
				"protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
				"protein_conformation__site_protein_conformation__site",
				Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
				annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role','ligand__properities__web_links__web_resource')),
				Prefetch("extra_proteins", queryset=StructureExtraProteins.objects.all().prefetch_related(
					'protein_conformation','wt_protein')),
				Prefetch("signprotcomplex_set", queryset=SignprotComplex.objects.all().prefetch_related('protein')))
		except Structure.DoesNotExist as e:
			pass

		structs_and_coverage = []
		for s in structures:
			structure_residues = Residue.objects.filter(protein_conformation=s.protein_conformation, protein_segment__isnull=False)
			coverage = round((len(structure_residues) / len(s.protein_conformation.protein.parent.sequence))*100)
			structs_and_coverage.append([s, coverage])
		context['structures'] = structs_and_coverage

		return context


class GProteinStructureBrowser(TemplateView):
	"""
	Fetching Structure data for browser
	"""
	template_name = "g_protein_structure_browser.html"

	def get_context_data (self, **kwargs):
		# Fetch g prot - receptor compleces
		context = super(GProteinStructureBrowser, self).get_context_data(**kwargs)
		complexstructs = SignprotComplex.objects.filter(protein__family__slug__startswith='100')
		try:
			context['structures'] = Structure.objects.filter(id__in=complexstructs.values_list('structure', flat=True)).select_related(
				"state",
				"pdb_code__web_resource",
				"protein_conformation__protein__species",
				"protein_conformation__protein__source",
				"protein_conformation__protein__family__parent__parent__parent",
				"publication__web_link__web_resource").prefetch_related(
				"stabilizing_agents", "construct__crystallization__crystal_method",
				"protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
				"protein_conformation__site_protein_conformation__site",
				Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
				annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role','ligand__properities__web_links__web_resource')),
				Prefetch("extra_proteins", queryset=StructureExtraProteins.objects.all().prefetch_related(
					'protein_conformation','wt_protein')),
				Prefetch("signprot_complex", queryset=SignprotComplex.objects.all().prefetch_related('protein')))
		except Structure.DoesNotExist as e:
			pass
		# Fetch non-complex g prot structures and filter for overlaps preferring SignprotComplex
		ncstructs = SignprotStructure.objects.filter(protein__family__slug__startswith='100').select_related(
				"protein__family",
				"pdb_code__web_resource",
				"publication__web_link__web_resource").prefetch_related(
				"stabilizing_agents",
				Prefetch("extra_proteins", queryset=SignprotStructureExtraProteins.objects.all().prefetch_related('wt_protein')))
		pdbs = []
		filtered_ncstructs = []
		for i in context['structures']:
			if i.pdb_code.index not in pdbs:
				pdbs.append(i.pdb_code.index)
		for i in ncstructs:
			if i.pdb_code.index not in pdbs:
				pdbs.append(i.pdb_code.index)
				filtered_ncstructs.append(i)
		context['structures'] = list(context['structures'])+list(filtered_ncstructs)
		return context


class ServeHomologyModels(TemplateView):

	template_name = "homology_models.html"
	def get_context_data(self, **kwargs):
		context = super(ServeHomologyModels, self).get_context_data(**kwargs)
		try:
			context['structure_model'] = StructureModel.objects.all().prefetch_related(
				"protein__family",
				"state",
				"protein__family__parent__parent__parent",
				"protein__species",
				"main_template__protein_conformation__protein__parent__family",
				"main_template__pdb_code")
		except StructureModel.DoesNotExist as e:
			pass

		return context


class ServeComplexModels(TemplateView):

	template_name = "complex_models.html"
	def get_context_data(self, **kwargs):
		context = super(ServeComplexModels, self).get_context_data(**kwargs)
		try:
			context['structure_complex_model'] = StructureComplexModel.objects.all().prefetch_related(
				"receptor_protein",
				"receptor_protein__family",
				"receptor_protein__family__parent__parent__parent",
				"receptor_protein__species",
				"sign_protein",
				"sign_protein__family",
				"sign_protein__family__parent__parent__parent",
				"main_template__protein_conformation__protein__parent__family",
				"main_template__pdb_code",
				"main_template__signprot_complex")
		except StructureComplexModel.DoesNotExist as e:
			pass

		return context


class ServeModelStatistics(TemplateView):

	template_name = "model_statistics.html"
	def get_context_data(self, **kwargs):
		context = super(ServeModelStatistics, self).get_context_data(**kwargs)
		smr = StructureModelRMSD.objects.all()
		try:
			context['structure_model_rmsds'] = smr.prefetch_related(
				"target_structure__protein_conformation__protein__parent",
				"target_structure__protein_conformation__protein__parent__family")
		except StructureModelRMSD.DoesNotExist as e:
			print(e)

		return context


def RedirectBrowser(request):
    response = redirect("/structure/")
    return response

def HomologyModelDetails(request, modelname, state):
	"""
	Show homology models details
	"""
	modelname = modelname

	color_palette = ["orange","cyan","yellow","lime","fuchsia","green","teal","olive","thistle","grey","chocolate","blue","red","pink","maroon"]

	model = StructureModel.objects.get(protein__entry_name=modelname, state__slug=state)
	model_main_template = model.main_template
	if model.protein.accession:
		residues = Residue.objects.filter(protein_conformation__protein=model.protein)
		a = Alignment()
		a.load_reference_protein(model.protein)
		a.load_proteins([model.main_template.protein_conformation.protein.parent])
		segs = ProteinSegment.objects.filter(id__in=residues.order_by("protein_segment__slug").distinct("protein_segment__slug").values_list("protein_segment", flat=True))
		a.load_segments(segs)
		a.build_alignment()
		a.calculate_similarity()
		main_template_seqsim = a.proteins[1].similarity
	else:
		residues = Residue.objects.filter(protein_conformation__protein=model.protein.parent)
		main_template_seqsim = 100
	rotamers = parse_model_statsfile(model.stats_text.stats_text, residues)
	version = model.version

	bb_temps, backbone_templates, r_temps, rotamer_templates, segments_out, bb_main, bb_alt, bb_none, sc_main, sc_alt, sc_none, template_list, colors = format_model_details(rotamers, model_main_template, color_palette)

	return render(request,'homology_models_details.html',{'model': model, 'modelname': modelname, 'rotamers': rotamers, 'backbone_templates': bb_temps, 'backbone_templates_number': len(backbone_templates),
														  'rotamer_templates': r_temps, 'rotamer_templates_number': len(rotamer_templates), 'color_residues': json.dumps(segments_out), 'bb_main': round(bb_main/len(rotamers)*100, 1),
														  'bb_alt': round(bb_alt/len(rotamers)*100, 1), 'bb_none': round(bb_none/len(rotamers)*100, 1), 'sc_main': round(sc_main/len(rotamers)*100, 1), 'sc_alt': round(sc_alt/len(rotamers)*100, 1),
														  'sc_none': round(sc_none/len(rotamers)*100, 1), 'main_template_seqsim': main_template_seqsim, 'template_list': template_list, 'model_main_template': model_main_template,
														  'state': state, 'version': version})

def RefinedModelDetails(request, pdbname):
	if len(pdbname) == 4:
		try:
			structure = Structure.objects.get(pdb_code__index=pdbname.upper())
			if structure.refined:
				complex_mod_details = SignprotComplex.objects.filter(structure=structure)
				if len(complex_mod_details) > 0:
					complex_mod = complex_mod_details.first()
					return ComplexModelDetails(request, pdbname.lower(), complex_mod.protein.entry_name)
				else:
					return HomologyModelDetails(request, pdbname.lower(), structure.state.slug)

			else:
				error = f"This structure ({pdbname}) does not have a refined model"
		except Structure.DoesNotExist:
			error = f"The structure {pdbname} does not exist in the GPCRdb"
	else:
		error = f"An incorrect PDB entry ({pdbname}) was specified"

	return HttpResponse(error)

def ComplexModelDetails(request, modelname, signprot):
	"""
	Show complex homology models details
	"""
	color_palette = ["orange","cyan","yellow","lime","fuchsia","limegreen","teal","olive","thistle","grey","chocolate","blue","red","pink","palegoldenrod","steelblue","tan","lightcoral","skyblue","papayawhip"]
	model = StructureComplexModel.objects.get(receptor_protein__entry_name=modelname, sign_protein__entry_name=signprot)
	main_template = model.main_template
	if model.receptor_protein.accession:
		receptor_residues = Residue.objects.filter(protein_conformation__protein=model.receptor_protein)
		signprot_residues = Residue.objects.filter(protein_conformation__protein=model.sign_protein)
		a = Alignment()
		a.load_reference_protein(model.receptor_protein)
		a.load_proteins([model.main_template.protein_conformation.protein.parent])
		segs = ProteinSegment.objects.filter(id__in=receptor_residues.order_by("protein_segment__slug").distinct("protein_segment__slug").values_list("protein_segment", flat=True))
		a.load_segments(segs)
		a.build_alignment()
		a.calculate_similarity()
		main_template_seqsim = a.proteins[1].similarity
	else:
		receptor_residues = Residue.objects.filter(protein_conformation__protein=model.receptor_protein.parent)
		signprot_residues = Residue.objects.filter(protein_conformation__protein=model.sign_protein)
		main_template_seqsim = 100
	receptor_rotamers, signprot_rotamers = parse_model_statsfile(model.stats_text.stats_text, receptor_residues, signprot_residues)

	loop_segments = ProteinSegment.objects.filter(category='loop', proteinfamily='Alpha')


	signprot_template = SignprotComplex.objects.get(structure=main_template).protein
	bb_temps, backbone_templates, r_temps, rotamer_templates, segments_out, bb_main, bb_alt, bb_none, sc_main, sc_alt, sc_none, template_list, colors = format_model_details(receptor_rotamers, main_template, color_palette, chain='R')
	signprot_color_palette = [i for i in color_palette if i not in list(colors.values())]

	bb_temps2, backbone_templates2, r_temps2, rotamer_templates2, segments_out2, bb_main2, bb_alt2, bb_none2, sc_main2, sc_alt2, sc_none2, template_list2, colors2 = format_model_details(signprot_rotamers, main_template, signprot_color_palette, chain='A', used_colors=colors)

	gp = GProteinAlignment()
	gp.run_alignment(model.sign_protein, signprot_template, calculate_similarity=True)

	for n in bb_temps2.values():
		for s in n:
			if s.protein_conformation.protein.parent not in bb_temps:
				bb_temps[s.protein_conformation.protein.parent] = [s]
			else:
				if s not in bb_temps[s.protein_conformation.protein.parent]:
					bb_temps[s.protein_conformation.protein.parent].append(s)
					break
	return render(request,'complex_models_details.html',{'model': model, 'modelname': modelname, 'signprot': signprot, 'signprot_template': signprot_template, 'receptor_rotamers': receptor_rotamers, 'signprot_rotamers': signprot_rotamers, 'backbone_templates': bb_temps, 'backbone_templates_number': len(backbone_templates),
														 'rotamer_templates': r_temps, 'rotamer_templates_number': len(rotamer_templates), 'color_residues': json.dumps(segments_out), 'bb_main': round(bb_main/len(receptor_rotamers)*100, 1),
														 'bb_alt': round(bb_alt/len(receptor_rotamers)*100, 1), 'bb_none': round(bb_none/len(receptor_rotamers)*100, 1), 'sc_main': round(sc_main/len(receptor_rotamers)*100, 1),
														 'sc_alt': round(sc_alt/len(receptor_rotamers)*100, 1), 'sc_none': round(sc_none/len(receptor_rotamers)*100, 1), 'main_template_seqsim': main_template_seqsim,
														 'template_list': template_list, 'model_main_template': main_template, 'state': None, 'signprot_sim': int(gp.proteins[1].similarity),
														 'signprot_color_residues': json.dumps(segments_out2), 'loop_segments': loop_segments})#, 'delta_distance': delta_distance})

def parse_model_statsfile(statstext, receptor_residues, signprot_residues=None):
	receptor_rotamers = []
	receptor_residues_dict = {r.sequence_number:r for r in receptor_residues}
	if signprot_residues:
		signprot_residues_dict = {r.sequence_number:r for r in signprot_residues}
		signprot_rotamers = []
		alpha_segments = ProteinSegment.objects.filter(proteinfamily='Alpha').values_list('slug', flat=True)
	structure_dict = {}

	for line in statstext.split('\n')[1:-1]:
		mr = ModelRotamer()
		split_line = line.split(',')
		del split_line[2]
		if (len(split_line) == 4):
			segment, seqnum, backbone_pdb, rotamer_pdb = split_line
		else:
			del split_line[2]
			segment, seqnum, backbone_pdb, rotamer_pdb = split_line

		if backbone_pdb not in structure_dict:
			if backbone_pdb=='None':
				backbone_struct = None
			else:
				backbone_struct = Structure.objects.get(pdb_code__index=backbone_pdb)
				structure_dict[backbone_pdb] = backbone_struct
		else:
			backbone_struct = structure_dict[backbone_pdb]
		if rotamer_pdb not in structure_dict:
			if rotamer_pdb=='None':
				rotamer_struct = None
			else:
				rotamer_struct = Structure.objects.get(pdb_code__index=rotamer_pdb)
				structure_dict[rotamer_pdb] = rotamer_struct
		else:
			rotamer_struct = structure_dict[rotamer_pdb]
		mr.backbone_template = backbone_struct
		mr.rotamer_template = rotamer_struct
		# Statsfile is for a complex model
		if signprot_residues and segment in alpha_segments:
			mr.residue = signprot_residues_dict[int(seqnum)]#signprot_residues.get(protein_segment__slug=segment, sequence_number=int(seqnum))
			signprot_rotamers.append(mr)
		elif segment in ['Beta','Gamma']:
			continue
		else:
			mr.residue = receptor_residues_dict[int(seqnum)]
			receptor_rotamers.append(mr)
	if signprot_residues:
		return receptor_rotamers, signprot_rotamers
	else:
		return receptor_rotamers

def format_model_details(rotamers, model_main_template, color_palette, chain=None, used_colors=None):
	backbone_templates, rotamer_templates = [],[]
	segments, segments_formatted, segments_out = {},{},{}
	bb_temps, r_temps = OrderedDict(), OrderedDict()
	bb_main, bb_alt, bb_none = 0,0,0
	sc_main, sc_alt, sc_none = 0,0,0

	for r in rotamers:
		if r.backbone_template not in backbone_templates and r.backbone_template!=None:
			backbone_templates.append(r.backbone_template)
			if r.backbone_template.protein_conformation.protein.parent not in bb_temps:
				bb_temps[r.backbone_template.protein_conformation.protein.parent] = [r.backbone_template]
			else:
				bb_temps[r.backbone_template.protein_conformation.protein.parent].append(r.backbone_template)
		if r.rotamer_template not in rotamer_templates and r.rotamer_template!=None:
			rotamer_templates.append(r.rotamer_template)
			if r.rotamer_template.protein_conformation.protein.parent not in r_temps:
				r_temps[r.rotamer_template.protein_conformation.protein.parent] = [r.rotamer_template]
			else:
				r_temps[r.rotamer_template.protein_conformation.protein.parent].append(r.rotamer_template)
		if r.backbone_template not in segments:
			segments[r.backbone_template] = [r.residue.sequence_number]
		else:
			segments[r.backbone_template].append(r.residue.sequence_number)
		if r.backbone_template==model_main_template:
			bb_main+=1
		elif r.backbone_template!=None:
			bb_alt+=1
		elif r.backbone_template==None:
			bb_none+=1
		if r.rotamer_template==model_main_template:
			sc_main+=1
		elif r.rotamer_template!=None:
			sc_alt+=1
		elif r.rotamer_template==None:
			sc_none+=1
	for s, nums in segments.items():
		for i, num in enumerate(nums):
			if i==0:
				segments_formatted[s] = [[num]]
			elif nums[i-1]!=num-1:
				if segments_formatted[s][-1][0]==nums[i-1]:
					segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1])
				else:
					segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1])
				segments_formatted[s].append([num])
				if i+1==len(segments[s]):
					segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], segments_formatted[s][-1][0])
			elif i+1==len(segments[s]):
				segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1]+1)
		if len(nums)==1:
			segments_formatted[s] = ['{}-{}'.format(segments_formatted[s][0][0], segments_formatted[s][0][0])]

	colors = OrderedDict([(model_main_template,"darkorchid"), (None,"white")])
	i = 0
	for s, nums in segments_formatted.items():
		if len(nums)>1:
			if chain:
				text = ':{} and ('.format(chain)
			else:
				text = ''
			for n in nums:
				text+='{} or '.format(n)
			if chain:
				segments_formatted[s] = text[:-4]+')'
			else:
				segments_formatted[s] = text[:-4]
		else:
			if chain:
				segments_formatted[s] = ':{} and ({})'.format(chain, segments_formatted[s][0])
			else:
				segments_formatted[s] = segments_formatted[s][0]
		if s==model_main_template:
			pass
		elif s==None:
			segments_out["white"] = segments_formatted[s]
		else:
			if used_colors:
				if s in used_colors:
					segments_out[used_colors[s]] = segments_formatted[s]
					colors[s] = used_colors[s]
				else:
					segments_out[color_palette[i]] = segments_formatted[s]
					colors[s] = color_palette[i]
			else:
				segments_out[color_palette[i]] = segments_formatted[s]
				colors[s] = color_palette[i]
		i+=1
	template_list = []
	for b, temps in bb_temps.items():
		for i, t in enumerate(temps):
			t.color = colors[t]
			bb_temps[b][i] = t
			template_list.append(t.pdb_code.index)

	segments_out = [[i,j] for i,j in segments_out.items()]

	return bb_temps, backbone_templates, r_temps, rotamer_templates, segments_out, bb_main, bb_alt, bb_none, sc_main, sc_alt, sc_none, template_list, colors

def ServeHomModDiagram(request, modelname, state):
	model=StructureModel.objects.filter(protein__entry_name=modelname, state__slug=state)
	if model.exists():
		model=model.get()
	else:
		 quit() #quit!

	if model.pdb_data is None:
		quit()

	response = HttpResponse(model.pdb_data.pdb, content_type='text/plain')
	return response

def ServeComplexModDiagram(request, modelname, signprot):
	model=StructureComplexModel.objects.filter(receptor_protein__entry_name=modelname, sign_protein__entry_name=signprot)
	if model.exists():
		model=model.get()
	else:
		 quit() #quit!

	if model.pdb_data is None:
		quit()

	response = HttpResponse(model.pdb_data.pdb, content_type='text/plain')
	return response

def StructureDetails(request, pdbname):
	"""
	Show structure details
	"""
	pdbname = pdbname.upper()
	structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name','structure_ligand_pair__pdb_reference','structure_ligand_pair__annotated').filter(structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
	resn_list = ''

	main_ligand = 'None'
	for structure in structures:
		if structure['structure_ligand_pair__annotated']:
			resn_list += ",\""+structure['structure_ligand_pair__pdb_reference']+"\""
			main_ligand = structure['structure_ligand_pair__pdb_reference']

	crystal = Structure.objects.get(pdb_code__index=pdbname)
	ligands = StructureLigandInteraction.objects.filter(structure=crystal, annotated=True)
	p = Protein.objects.get(protein=crystal.protein_conformation.protein)
	residues = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number')

	# positioning data
	sv = StructureVectors.objects.filter(structure=crystal)
	translation = center_axis = ""
	if sv.exists():
		sv = sv.get()
		translation = sv.translation
		center_axis = sv.center_axis

	# Check if the structure is in complex with a signaling protein
	signaling_complex = SignprotComplex.objects.filter(structure=crystal).count() > 0

	# GN list
	only_gns = list(crystal.protein_conformation.residue_set.exclude(generic_number=None).values_list('protein_segment__slug','sequence_number','generic_number__label','display_generic_number__label').all())
	gn_list = [x[1] for x in only_gns]
	filter_tm1 = [x[1] for x in only_gns if x[2] == "1x46"]
	ref_tm1 = ""
	if len(filter_tm1) > 0:
		ref_tm1 = filter_tm1[0]

	return render(request,'structure_details.html',{'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein':p, 'residues':residues, 'annotated_resn': resn_list, 'main_ligand': main_ligand, 'ligands': ligands, 'translation': translation, 'center_axis': center_axis, 'gn_list': gn_list, 'ref_tm1': ref_tm1, 'signaling_complex': signaling_complex})

def ServePdbDiagram(request, pdbname):
	structure=Structure.objects.filter(pdb_code__index=pdbname.upper())
	if structure.exists():
		structure=structure.get()
	else:
		 quit() #quit!

	if structure.pdb_data is None:
		quit()

	response = HttpResponse(structure.pdb_data.pdb, content_type='text/plain')
	return response


def ServeUprightPdbDiagram(request, pdbname):
	structure = Structure.objects.filter(pdb_code__index=pdbname.upper())
	if structure.exists():
		structure = structure.get()
	else:
		 quit() #quit!

	if structure.pdb_data is None:
		quit()

	sv = StructureVectors.objects.filter(structure=structure)
	struct = translation = center_axis = ""
	if sv.exists():
		sv = sv.get()
		translation = json.loads(sv.translation)
		center_axis = json.loads(sv.center_axis)

	# Load structure
	parser = PDBParser(QUIET=True)
	with io.StringIO(structure.pdb_data.pdb) as f:
		struct = parser.get_structure(structure.pdb_code.index, f)

	# Neutral references
	translation_neutral = np.array((0, 0, 0), 'f')
	qn_neutral = rotmat(Vector(0, 1, 0), Vector(0, 1, 0))

	# Calculate translation
	translation_matrix = np.array((translation[0], translation[1], translation[2]), 'f')

	# Calculate rotation
	v1 = Vector(0, 1, 0)
	v2 = Vector(center_axis[0], center_axis[1], center_axis[2])
	v3 = Vector(-1, 0, 0)
	qn_upright = rotmat(v1, v2)

	# Calculate rotation around axis (TM1 placement)
	ref_tm1 = Residue.objects.filter(protein_conformation=structure.protein_conformation, generic_number_id__label="1x46")
	if ref_tm1.count() > 0:
		ref_tm1 = ref_tm1.first().sequence_number
		tm1_coord = struct[0][structure.preferred_chain][ref_tm1]["CA"].get_vector()
		ref_vector = Vector(tm1_coord[0], 0, tm1_coord[2]) #height position doesn't matter
		ref_vector = ref_vector.normalized()
		qn_axis = rotmat(ref_vector, v3)


	# Apply transformations
	for atom in struct.get_atoms():
		# First apply translation before rotation
		atom.transform(qn_neutral, translation_matrix)
		# Second rotate - calculate rotation and apply
		atom.transform(qn_upright, translation_neutral)
		# Finally, rotate around axis to place helices correctly
		atom.transform(qn_axis, translation_neutral)

	# Save PDB transformed structure
	out_stream = StringIO()
	pdb_out = PDBIO()
	pdb_out.set_structure(struct)
	pdb_out.save(out_stream)

	# Send as response
	return HttpResponse(out_stream.getvalue(), content_type='chemical/x-pdb')

def ServePdbLigandDiagram(request,pdbname,ligand):
	pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(Q(ligand__properities__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
	response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
	return response

class StructureStatistics(TemplateView):
	"""
	So not ready that EA wanted to publish it.
	"""
	template_name = 'structure_statistics.html'
	origin = 'structure'

	def get_context_data (self, **kwargs):
		context = super().get_context_data(**kwargs)
		families = ProteinFamily.objects.all()
		lookup = {}
		for f in families:
			lookup[f.slug] = f.name

		all_structs = Structure.objects.all().prefetch_related('protein_conformation__protein__family')
		noncomplex_gprots = SignprotStructure.objects.filter(protein__family__slug__startswith='100')
		all_complexes = all_structs.exclude(ligands=None)
		all_gprots = all_structs.filter(id__in=SignprotComplex.objects.filter(protein__family__slug__startswith='100').values_list("structure__id", flat=True))
		###### these are query sets for G-Prot Structure Statistics
		all_g_A_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='001')
		all_g_B1_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='002')
		all_g_B2_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='003')
		all_g_C_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='004')
		all_g_D1_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='005')
		all_g_F_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='006')
		all_g_T2_complexes = all_gprots.filter(protein_conformation__protein__family__slug__startswith='007')
		######
		all_active = all_structs.filter(protein_conformation__state__slug = 'active')

		years = self.get_years_range(list(set([x.publication_date.year for x in all_structs])))

		unique_structs = Structure.objects.order_by('protein_conformation__protein__family__name', 'state',
			'publication_date', 'resolution').distinct('protein_conformation__protein__family__name').prefetch_related('protein_conformation__protein__family')
		unique_complexes = StructureLigandInteraction.objects.filter(annotated=True).distinct('ligand', 'structure__protein_conformation__protein__family')
		unique_gprots = unique_structs.filter(id__in=SignprotComplex.objects.filter(protein__family__slug__startswith='100').values_list("structure__id", flat=True))
		# unique_g_A_complexes = unique_gprots.filter(protein_conformation__protein__family__slug__startswith='001')
		unique_g_A_complexes = all_g_A_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_B1_complexes = all_g_B1_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_B2_complexes = all_g_B2_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_C_complexes = all_g_C_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_D1_complexes = all_g_D1_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_F_complexes = all_g_F_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_g_T2_complexes = all_g_T2_complexes.annotate(distinct_name=Concat('signprot_complex__protein__family__parent__name', 'protein_conformation__protein__family__name', output_field=TextField())).order_by('distinct_name').distinct('distinct_name')
		unique_active = unique_structs.filter(protein_conformation__state__slug = 'active')
		tree_dots_data = {}
		tree_dots_data = self.grab_matches(unique_g_A_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_B1_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_B2_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_C_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_D1_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_F_complexes, tree_dots_data)
		tree_dots_data = self.grab_matches(unique_g_T2_complexes, tree_dots_data)
		#Stats
		struct_count = Structure.objects.all().annotate(Count('id'))
		struct_lig_count = Structure.objects.exclude(ligands=None)

		context['all_structures'] = len(all_structs)
		context['all_structures_by_class'] = self.count_by_class(all_structs, lookup)
		context['all_complexes'] = len(all_complexes)
		context['all_complexes_by_class'] = self.count_by_class(all_complexes, lookup)
		context['all_gprots'] = len(all_gprots)
		context['all_gprots_by_class'] = self.count_by_class(all_gprots, lookup)
		context['all_gprots_by_gclass'] = self.count_by_gclass(all_gprots, lookup)
		context['noncomplex_gprots_by_gclass'] = self.count_by_gclass(noncomplex_gprots, lookup, True)
		context['noncomplex_gprots'] = len(noncomplex_gprots)

		context['gA_complexes'] = zip(list(self.count_by_gclass(all_g_A_complexes, lookup).items()), list(self.count_by_gclass(unique_g_A_complexes, lookup).items()))
		context['all_g_A_complexes'] = len(all_g_A_complexes)
		context['unique_g_A_complexes'] = len(unique_g_A_complexes)
		context['gB1_complexes'] = zip(list(self.count_by_gclass(all_g_B1_complexes, lookup).items()), list(self.count_by_gclass(unique_g_B1_complexes, lookup).items()))
		context['all_g_B1_complexes'] = len(all_g_B1_complexes)
		context['unique_g_B1_complexes'] = len(unique_g_B1_complexes)
		context['gB2_complexes'] = zip(list(self.count_by_gclass(all_g_B2_complexes, lookup).items()), list(self.count_by_gclass(unique_g_B2_complexes, lookup).items()))
		context['all_g_B2_complexes'] = len(all_g_B2_complexes)
		context['unique_g_B2_complexes'] = len(unique_g_B2_complexes)
		context['gC_complexes'] = zip(list(self.count_by_gclass(all_g_C_complexes, lookup).items()), list(self.count_by_gclass(unique_g_C_complexes, lookup).items()))
		context['all_g_C_complexes'] = len(all_g_C_complexes)
		context['unique_g_C_complexes'] = len(unique_g_C_complexes)
		context['gD1_complexes'] = zip(list(self.count_by_gclass(all_g_D1_complexes, lookup).items()), list(self.count_by_gclass(unique_g_D1_complexes, lookup).items()))
		context['all_g_D1_complexes'] = len(all_g_D1_complexes)
		context['unique_g_D1_complexes'] = len(unique_g_D1_complexes)
		context['gF_complexes'] = zip(list(self.count_by_gclass(all_g_F_complexes, lookup).items()), list(self.count_by_gclass(unique_g_F_complexes, lookup).items()))
		context['all_g_F_complexes'] = len(all_g_F_complexes)
		context['unique_g_F_complexes'] = len(unique_g_F_complexes)
		context['gT2_complexes'] = zip(list(self.count_by_gclass(all_g_T2_complexes, lookup).items()), list(self.count_by_gclass(unique_g_T2_complexes, lookup).items()))
		context['all_g_T2_complexes'] = len(all_g_T2_complexes)
		context['unique_g_T2_complexes'] = len(unique_g_T2_complexes)

		context['all_active'] = len(all_active)
		context['all_active_by_class'] = self.count_by_class(all_active, lookup)
		context['unique_structures'] = len(unique_structs)
		context['unique_structures_by_class'] = self.count_by_class(unique_structs, lookup)
		context['unique_complexes'] = len(unique_complexes)
		context['unique_complexes_by_class'] = self.count_by_class([x.structure for x in unique_complexes], lookup)
		context['unique_gprots'] = len(unique_gprots)
		context['unique_gprots_by_gclass'] = self.count_by_gclass(unique_gprots, lookup)
		context['unique_gprots_by_class'] = self.count_by_class(unique_gprots, lookup)
		context['unique_active'] = len(unique_active)
		context['unique_active_by_class'] = self.count_by_class(unique_active, lookup)
		context['release_notes'] = ReleaseNotes.objects.all()[0]
		context['latest_structure'] = Structure.objects.latest('publication_date').publication_date

		context['chartdata'] = self.get_per_family_cumulative_data_series(years, unique_structs, lookup)
		context['chartdata_y'] = self.get_per_family_data_series(years, unique_structs, lookup)
		context['chartdata_all'] = self.get_per_family_cumulative_data_series(years, all_structs, lookup)
		context['chartdata_reso'] = self.get_resolution_coverage_data_series(all_structs)

		context['chartdata_class'] = self.get_per_class_cumulative_data_series(years, unique_structs, lookup)
		context['chartdata_class_y'] = self.get_per_class_data_series(years, unique_structs, lookup)
		context['chartdata_class_all'] = self.get_per_class_cumulative_data_series(years, all_structs, lookup)

		context['tree_dots_data'] = tree_dots_data
		#context['coverage'] = self.get_diagram_coverage()
		#{
		#    'depth': 3,
		#    'anchor': '#crystals'}
		# relabeling table columns for sake of consistency
		for key in list(context['unique_structures_by_class'].keys()):
			context['unique_structures_by_class'][key.replace('Class','')] = context['unique_structures_by_class'].pop(key)
		for key in list(context['all_structures_by_class'].keys()):
			context['all_structures_by_class'][key.replace('Class','')] = context['all_structures_by_class'].pop(key)
		context['total_gprots_by_gclass'] = []
		for key in context['all_gprots_by_gclass']:
			context['total_gprots_by_gclass'].append(context['all_gprots_by_gclass'][key] + context['noncomplex_gprots_by_gclass'][key])
		context['total_gprots'] = sum(context['total_gprots_by_gclass'])

		tree = PhylogeneticTreeGenerator()
		class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
		context['class_a_options'] = deepcopy(tree.d3_options)
		context['class_a_options']['anchor'] = 'class_a'
		context['class_a_options']['leaf_offset'] = 50
		context['class_a_options']['label_free'] = []
        # section to remove Orphan from Class A tree and apply to a different tree
		whole_class_a = class_a_data.get_nodes_dict('crystals')
		for item in whole_class_a['children']:
			if item['name'] == 'Orphan':
				orphan_data = OrderedDict([('name', ''), ('value', 3000), ('color', ''), ('children',[item])])
				whole_class_a['children'].remove(item)
				break
		context['class_a'] = json.dumps(whole_class_a)
		class_b1_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
		context['class_b1_options'] = deepcopy(tree.d3_options)
		context['class_b1_options']['anchor'] = 'class_b1'
		context['class_b1_options']['branch_trunc'] = 60
		context['class_b1_options']['label_free'] = [1,]
		context['class_b1'] = json.dumps(class_b1_data.get_nodes_dict('crystals'))
		class_b2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
		context['class_b2_options'] = deepcopy(tree.d3_options)
		context['class_b2_options']['anchor'] = 'class_b2'
		context['class_b2_options']['label_free'] = [1,]
		context['class_b2'] = json.dumps(class_b2_data.get_nodes_dict('crystals'))
		class_c_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
		context['class_c_options'] = deepcopy(tree.d3_options)
		context['class_c_options']['anchor'] = 'class_c'
		context['class_c_options']['branch_trunc'] = 50
		context['class_c_options']['label_free'] = [1,]
		context['class_c'] = json.dumps(class_c_data.get_nodes_dict('crystals'))
		class_f_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
		context['class_f_options'] = deepcopy(tree.d3_options)
		context['class_f_options']['anchor'] = 'class_f'
		context['class_f_options']['label_free'] = [1,]
		#json.dump(class_f_data.get_nodes_dict('crystalized'), open('tree_test.json', 'w'), indent=4)
		context['class_f'] = json.dumps(class_f_data.get_nodes_dict('crystals'))
		class_t2_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class T (Taste 2)'))
		context['class_t2_options'] = deepcopy(tree.d3_options)
		context['class_t2_options']['anchor'] = 'class_t2'
		context['class_t2_options']['label_free'] = [1,]
		context['class_t2'] = json.dumps(class_t2_data.get_nodes_dict('crystals'))
		# definition of the class a orphan tree
		context['orphan_options'] = deepcopy(tree.d3_options)
		context['orphan_options']['anchor'] = 'orphan'
		context['orphan_options']['label_free'] = [1,]
		context['orphan'] = json.dumps(orphan_data)
		whole_receptors = Protein.objects.prefetch_related("family", "family__parent__parent__parent").filter(sequence_type__slug="wt", family__slug__startswith="00")
		whole_rec_dict = {}
		for rec in whole_receptors:
			rec_uniprot = rec.entry_short()
			rec_iuphar = rec.family.name.replace("receptor", '').replace("<i>","").replace("</i>","").strip()
			if (rec_iuphar[0].isupper()) or (rec_iuphar[0].isdigit()):
				whole_rec_dict[rec_uniprot] = [rec_iuphar]
			else:
				whole_rec_dict[rec_uniprot] = [rec_iuphar.capitalize()]

		context["whole_receptors"] = json.dumps(whole_rec_dict)
		context["page"] = self.origin
		return context

	def get_families_dict(self, queryset, lookup):

		families = []
		for s in queryset:
			fid = s.protein_conformation.protein.family.slug.split("_")
			fname = lookup[fid[0]+"_"+fid[1]]
			cname = lookup[fid[0]]
			if fname not in families:
				families.append(fname)
		return families

	def grab_matches(self, queryset, output):

		#Grab data from queryset
		for s in queryset:
			gprot = s.signprot_complex.protein.family.parent.name
			receptor = s.protein_conformation.protein.parent.entry_short()
			if receptor in output.keys():
				output[receptor].append(gprot)
			else:
				output[receptor] = []
				output[receptor].append(gprot)

		return output

	def count_by_class(self, queryset, lookup):

		#Ugly walkaround
		classes = [lookup[x] for x in reversed(['001', '002', '003', '004', '005', '006', '007'])]
		records = []
		for s in queryset:
			fid = s.protein_conformation.protein.family.slug.split("_")
			cname = lookup[fid[0]]
			records.append(cname)

		tmp = OrderedDict()
		for x in sorted(classes):
			tmp[x] = records.count(x)

		return tmp

	def count_by_gclass(self, queryset, lookup, nc=False):

		#Ugly walkaround
		classes = [lookup[x] for x in ['100_001_001', '100_001_002', '100_001_003', '100_001_004', '100_001_005']]
		translate = {'Gs':'G<sub>s</sub>', 'Gi/o':'G<sub>i/o</sub>', 'Gq/11':'G<sub>q/11</sub>', 'G12/13':'G<sub>12/13</sub>', 'GPa1 family':'GPa1'}
		records = []
		if nc == False:
			for s in queryset:
				fid = s.signprot_complex.protein.family.parent.slug
				cname = lookup[fid]
				records.append(translate[cname])
		else:
			for s in queryset:
				fid = s.protein.family.parent.slug
				cname = lookup[fid]
				records.append(translate[cname])

		tmp = OrderedDict()
		for x in classes:
			tmp[translate[x]] = records.count(translate[x])

		return tmp

	def get_years_range(self, years_list):

		min_y = min(years_list)
		max_y = max(years_list)
		return range(min_y, max_y+1)

	def get_per_class_data_series(self, years, structures, lookup):
		"""
		Prepare data for multiBarGraph of unique crystallized receptors grouped by class. Returns data series for django-nvd3 wrapper.
		"""
		classes = [lookup[x] for x in ['001', '002', '003', '004', '005', '006', '007']]
		series = []
		data = {}
		for year in years:
			for prot_class in classes:
				if prot_class not in data.keys():
					data[prot_class] = []
				count = 0
				for structure in structures:
					fid = structure.protein_conformation.protein.family.slug.split("_")
					# if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
					if lookup[fid[0]] == prot_class and structure.publication_date.year == year:
						count += 1
				data[prot_class].append(count)
		for prot_class in classes:
			series.append({"values":
				[OrderedDict({
					'x': years[i],
					'y': j
					}) for i, j in enumerate(data[prot_class])],
				"key": prot_class,
				"yAxis": "1"})
		return json.dumps(series)

	def get_per_family_data_series(self, years, structures, lookup):
		"""
		Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
		"""
		families = self.get_families_dict(structures, lookup)
		series = []
		data = {}
		for year in years:
			for family in families:
				if family not in data.keys():
					data[family] = []
				count = 0
				for structure in structures:
					fid = structure.protein_conformation.protein.family.slug.split("_")
					# if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
					if lookup[fid[0]+"_"+fid[1]] == family and structure.publication_date.year == year:
						count += 1
				data[family].append(count)
		for family in families:
			series.append({"values":
				[OrderedDict({
					'x': years[i],
					'y': j
					}) for i, j in enumerate(data[family])],
				"key": family,
				"yAxis": "1"})
		return json.dumps(series)

	def get_per_class_cumulative_data_series(self, years, structures, lookup):
		"""
		Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
		"""
		classes =  [lookup[x] for x in ['001', '002', '003', '004', '005', '006', '007']]
		series = []
		data = {}
		for year in years:
			for prot_class in classes:
				if prot_class not in data.keys():
					data[prot_class] = []
				count = 0
				for structure in structures:
					fid = structure.protein_conformation.protein.family.slug.split("_")
					# if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
					if lookup[fid[0]] == prot_class and structure.publication_date.year == year:
						count += 1
				if len(data[prot_class]) > 0:
					data[prot_class].append(count + data[prot_class][-1])
				else:
					data[prot_class].append(count)
		for prot_class in classes:
			series.append({"values":
				[OrderedDict({
					'x': years[i],
					'y': j
					}) for i, j in enumerate(data[prot_class])],
				"key": prot_class,
				"yAxis": "1"})
		return json.dumps(series)


	def get_per_family_cumulative_data_series(self, years, structures, lookup):
		"""
		Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
		"""
		families = self.get_families_dict(structures, lookup)
		series = []
		data = {}
		for year in years:
			for family in families:
				if family not in data.keys():
					data[family] = []
				count = 0
				for structure in structures:
					fid = structure.protein_conformation.protein.family.slug.split("_")
					# if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
					if lookup[fid[0]+"_"+fid[1]] == family and structure.publication_date.year == year:
						count += 1
				if len(data[family]) > 0:
					data[family].append(count + data[family][-1])
				else:
					data[family].append(count)
		for family in families:
			series.append({"values":
				[{
					'x': years[i],
					'y': j
					} for i, j in enumerate(data[family])],
				"key": family,
				"yAxis": "1"})
		return json.dumps(series)


	def get_resolution_coverage_data_series(self, structures):
		"""
		Prepare data for multiBarGraph of resolution coverage of available crystal structures.
		"""
		#Resolutions boundaries
		reso_min = float(min([round(x.resolution, 1) for x in structures]))
		reso_max = float(max([round(x.resolution, 1) for x in structures]))
		step = (reso_max - reso_min)/10

		brackets = [reso_min + step*x for x in range(10)] + [reso_max]

		reso_count = []
		bracket_labels = []
		for idx, bracket in enumerate(brackets):
			if idx == 0:
				reso_count.append(len([x for x in structures if x.resolution <= bracket]))
				bracket_labels.append('< {:.1f}'.format(bracket))
			else:
				reso_count.append(len([x for x in structures if bracket-step < x.resolution <= bracket]))
				bracket_labels.append('{:.1f}-{:.1f}'.format(brackets[idx-1],bracket))

		return json.dumps([{"values": [{
					'x': bracket_labels[i],
					'y': j
					} for i, j in enumerate(reso_count)],
				"key": 'Resolution coverage',
				"yAxis": "1"}])

	def get_diagram_coverage(self):
		"""
		Prepare data for coverage diagram.
		"""

		families = ProteinFamily.objects.all()
		lookup = {}
		for f in families:
			lookup[f.slug] = f.name.replace("receptors","")

		class_proteins = Protein.objects.filter(family__slug__startswith="00", source__name='SWISSPROT').prefetch_related('family').order_by('family__slug')

		coverage = OrderedDict()

		temp = OrderedDict([
							('name',''),
							('interactions', 0),
							('receptor_i', 0) ,
							('mutations' , 0),
							('receptor_m', 0),
							('mutations_an' , 0),
							('receptor_m_an', 0),
							('receptor_t',0),
							('children', OrderedDict()) ,
							('fraction_i',0),
							('fraction_m',0),
							('fraction_m_an',0)
							])

		for p in class_proteins:
			fid = p.family.slug.split("_")
			if fid[0] not in coverage:
				coverage[fid[0]] = deepcopy(temp)
				coverage[fid[0]]['name'] = lookup[fid[0]]
			if fid[1] not in coverage[fid[0]]['children']:
				coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
				coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
			if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
			if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
				coverage[fid[0]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1


		class_interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__annotated=True).prefetch_related(
			'rotamer__residue__display_generic_number','interaction_type',
			'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
			'structure_ligand_pair__ligand__properities',
			)


		score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0} , 'interaction' : {},'mutation': {}}

		# Replace above as fractions etc is not required and it was missing xtals that didnt have interactions.
		unique_structs = list(Structure.objects.order_by('protein_conformation__protein__parent', 'state',
			'publication_date', 'resolution').distinct('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__family'))

		for s in unique_structs:
			fid = s.protein_conformation.protein.family.slug.split("_")
			coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_i'] = 1
			coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions'] += 1

		CSS_COLOR_NAMES = ["SteelBlue","SlateBlue","LightCoral","Orange","LightGreen","LightGray","PeachPuff","PaleGoldenRod"]

		tree = OrderedDict({'name':'GPCRs','children':[]})
		i = 0
		n = 0
		for c_v in coverage.values():
			c_v['name'] = c_v['name'].split("(")[0]
			if c_v['name'].strip() == 'Other GPCRs':
				continue
			children = []
			for lt_v in c_v['children'].values():
				if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
					continue
				children_rf = []
				for rf_v in lt_v['children'].values():
					rf_v['name'] = rf_v['name'].split("<")[0]
					if rf_v['name'].strip() == 'Class T (Taste 2)':
						continue
					children_r = []
					for r_v in rf_v['children'].values():
						r_v['color'] = CSS_COLOR_NAMES[i]
						r_v['sort'] = n
						children_r.append(r_v)
						n += 1
					rf_v['children'] = children_r
					rf_v['sort'] = n
					rf_v['color'] = CSS_COLOR_NAMES[i]
					children_rf.append(rf_v)
				lt_v['children'] = children_rf
				lt_v['sort'] = n
				lt_v['color'] = CSS_COLOR_NAMES[i]
				children.append(lt_v)
			c_v['children'] = children
			c_v['sort'] = n
			c_v['color'] = CSS_COLOR_NAMES[i]
			tree['children'].append(c_v)
			i += 1

		return json.dumps(tree)

	@staticmethod
	def get_diagram_crystals():
		"""
		Prepare data for coverage diagram.
		"""

		crystal_proteins = [x.protein_conformation.protein.parent for x in Structure.objects.order_by('protein_conformation__protein__parent', 'state',
			'publication_date', 'resolution').distinct('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__parent__family')]

		families = []
		for cryst_prot in crystal_proteins:
			families.append(cryst_prot.family)
			tmp = cryst_prot.family
			while tmp.parent is not None:
				tmp = tmp.parent
				families.append(tmp)
		lookup = {}
		for f in families:
			lookup[f.slug] = f.name.replace("receptors","")

		coverage = OrderedDict()
		temp = OrderedDict([
							('name',''),
							('interactions', 0),
							('receptor_i', 0) ,
							('mutations' , 0),
							('receptor_m', 0),
							('mutations_an' , 0),
							('receptor_m_an', 0),
							('receptor_t',0),
							('children', OrderedDict()) ,
							('fraction_i',0),
							('fraction_m',0),
							('fraction_m_an',0)
							])

		for p in crystal_proteins:
			fid = p.family.slug.split("_")
			if fid[0] not in coverage:
				coverage[fid[0]] = deepcopy(temp)
				coverage[fid[0]]['name'] = lookup[fid[0]]
			if fid[1] not in coverage[fid[0]]['children']:
				coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
				coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
			if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
			if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
				coverage[fid[0]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
				coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1


		class_interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__annotated=True).prefetch_related(
			'rotamer__residue__display_generic_number','interaction_type',
			'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
			'structure_ligand_pair__ligand__properities',
			)

		score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0} , 'interaction' : {},'mutation': {}}

		# Replace above as fractions etc is not required and it was missing xtals that didnt have interactions.
		unique_structs = list(Structure.objects.order_by('protein_conformation__protein__family__name', 'state',
			'publication_date', 'resolution').distinct('protein_conformation__protein__family__name').prefetch_related('protein_conformation__protein__family'))

		for s in unique_structs:
			fid = s.protein_conformation.protein.family.slug.split("_")
			coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_i'] = 1
			coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions'] += 1

		CSS_COLOR_NAMES = ["SteelBlue","SlateBlue","LightCoral","Orange","LightGreen","LightGray","PeachPuff","PaleGoldenRod"]

		tree = OrderedDict({'name':'GPCRs','children':[]})
		i = 0
		n = 0
		for c,c_v in coverage.items():
			c_v['name'] = c_v['name'].split("(")[0]
			if c_v['name'].strip() == 'Other GPCRs':
				continue
			children = []
			for lt,lt_v in c_v['children'].items():
				if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
					continue
				children_rf = []
				for rf,rf_v in lt_v['children'].items():
					rf_v['name'] = rf_v['name'].split("<")[0]
					if rf_v['name'].strip() == 'Class T (Taste 2)':
						continue
					children_r = []
					for r,r_v in rf_v['children'].items():
						r_v['color'] = CSS_COLOR_NAMES[i]
						r_v['sort'] = n
						children_r.append(r_v)
						n += 1
					rf_v['children'] = children_r
					rf_v['sort'] = n
					rf_v['color'] = CSS_COLOR_NAMES[i]
					children_rf.append(rf_v)
				lt_v['children'] = children_rf
				lt_v['sort'] = n
				lt_v['color'] = CSS_COLOR_NAMES[i]
				children.append(lt_v)
			c_v['children'] = children
			c_v['sort'] = n
			c_v['color'] = CSS_COLOR_NAMES[i]
			tree['children'].append(c_v)
			i += 1

		return json.dumps(tree)


class GenericNumberingIndex(TemplateView):
	"""
	Starting page of generic numbering assignment workflow.
	"""
	template_name = 'common_structural_tools.html'

	#Left panel
	step = 1
	number_of_steps = 2
	documentation_url = settings.DOCUMENTATION_URL
	docs = 'structures.html#pdb-file-residue-numbering'
	title = "UPLOAD A PDB FILE"
	description = """
	Upload a pdb file to be annotated with generic numbers from GPCRdb.

	The numbers can be visualized in molecular viewers such as PyMOL, with scripts available with the output files.

	Once you have selected all your targets, click the green button.
	"""

	#Input file form data
	header = "Select a file to upload:"
	upload_form_data = {
		"pdb_file": forms.FileField(),
		}
	form_code = forms.Form()
	form_code.fields = upload_form_data
	form_id = 'gn_pdb_file'
	url = '/structure/generic_numbering_results'
	mid_section = "upload_file_form.html"
	form_height = 200
	#Buttons
	buttons = {
		'continue' : {
			'label' : 'Assign generic numbers',
			'color' : 'success',
			},
		}


	def get_context_data (self, **kwargs):

		context = super(GenericNumberingIndex, self).get_context_data(**kwargs)
		# get attributes of this class and add them to the context
		context['form_code'] = str(self.form_code)
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		return context


#Class rendering results from generic numbers assignment
class GenericNumberingResults(TemplateView):

	template_name = 'common_structural_tools.html'

	#Left panel
	step = 1
	number_of_steps = 2
	title = "SELECT SUBSTRUCTURE"
	description = 'Download the desired substructures.'
	#Mid section
	mid_section = 'gn_results.html'
	#Buttons - none


	def post (self, request, *args, **kwargs):

		generic_numbering = GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore")))
		out_struct = generic_numbering.assign_generic_numbers()
		out_stream = StringIO()
		io = PDBIO()
		io.set_structure(out_struct)
		io.save(out_stream)
		if len(out_stream.getvalue()) > 0:
			request.session['gn_outfile'] = out_stream
			request.session['gn_outfname'] = request.FILES['pdb_file'].name
			self.success = True
		else:
			self.input_file = request.FILES['pdb_file'].name
			self.success = False

		context = super(GenericNumberingResults, self).get_context_data(**kwargs)
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		return render(request, self.template_name, context)



class GenericNumberingSelection(AbsSegmentSelection):
	"""
	Segment selection for download of annotated substructure.
	"""

	step = 2
	number_of_steps = 2

	docs = 'structures.html#pdb-file-residue-numbering'

	#Mid section
	#mid_section = 'segment_selection.html'

	#Right panel
	segment_list = True
	buttons = {
		'continue': {
			'label': 'Download substructure',
			'url': '/structure/generic_numbering_results/substr',
			'color': 'success',
		},
	}
	# OrderedDict to preserve the order of the boxes
	selection_boxes = OrderedDict([('reference', False),
		('targets', False),
		('segments', True),])


	def get_context_data(self, **kwargs):

		context = super(GenericNumberingSelection, self).get_context_data(**kwargs)

		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		context['selection'] = {}
		context['selection']['site_residue_groups'] = selection.site_residue_groups
		context['selection']['active_site_residue_group'] = selection.active_site_residue_group
		for selection_box, include in self.selection_boxes.items():
			if include:
				context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

		# get attributes of this class and add them to the context
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		return context


class GenericNumberingDownload(View):
	"""
	Serve the (sub)structure depending on user's choice.
	"""
	def get(self, request, *args, **kwargs):

		if self.kwargs['substructure'] == 'custom':
			return HttpResponseRedirect('/structure/generic_numbering_selection')

		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)
		out_stream = StringIO()
		io = PDBIO()
		request.session['gn_outfile'].seek(0)
		gn_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(request.session['gn_outfname'], request.session['gn_outfile'])[0]

		if self.kwargs['substructure'] == 'full':
			io.set_structure(gn_struct)
			io.save(out_stream)

		if self.kwargs['substructure'] == 'substr':
			io.set_structure(gn_struct)
			io.save(out_stream, GenericNumbersSelector(parsed_selection=SelectionParser(selection)))

		root, ext = os.path.splitext(request.session['gn_outfname'])
		response = HttpResponse(content_type="chemical/x-pdb")
		response['Content-Disposition'] = 'attachment; filename="{}_GPCRDB.pdb"'.format(root)
		response.write(out_stream.getvalue())

		return response


#==============================================================================

#========================Superposition of structures===========================
#Class for starting page of superposition workflow
class SuperpositionWorkflowIndex(TemplateView):

	template_name = "common_structural_tools.html"

	#Left panel
	step = 1
	number_of_steps = 3
	documentation_url = settings.DOCUMENTATION_URL
	docs = 'structures.html#structure-superposition'
	title = "UPLOAD YOUR FILES"
	description = """
	Upload a pdb file for reference structure, and one or more files that will be superposed. You can also select the structures from crystal structure browser.

	Once you have uploaded/selected all your targets, click the green button.
	"""

	header = "Upload or select your structures:"
	#
	upload_form_data = OrderedDict([
		('ref_file', forms.FileField(label="Reference structure")),
		('alt_files', MultiFileField(label="Structure(s) to superpose", max_num=10, min_num=1)),
		#('exclusive', forms.BooleanField(label='Download only superposed subset of atoms', widget=forms.CheckboxInput())),
		])
	form_code = forms.Form()
	form_code.fields = upload_form_data
	form_id = 'superpose_files'
	url = '/structure/superposition_workflow_selection'
	mid_section = 'superposition_workflow_upload_file_form.html'

	#Buttons
	buttons = {
		'continue' : {
			'label' : 'Select segments',
			'color' : 'success',
			}
		}

	# OrderedDict to preserve the order of the boxes
	selection_boxes = OrderedDict([('reference', True),
		('targets', True),
		('segments', False)])

	def get_context_data (self, **kwargs):
		context = super(SuperpositionWorkflowIndex, self).get_context_data(**kwargs)

		# get selection from session and add to context
		# get simple selection from session
		simple_selection = self.request.session.get('selection', False)
		# print(simple_selection)
		# create full selection and import simple selection (if it exists)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)
		# print(self.kwargs.keys())
		#Clearing selections for fresh run
		if 'clear' in self.kwargs.keys():
			selection.clear('reference')
			selection.clear('targets')
			selection.clear('segments')
			if 'alt_files' in self.request.session.keys():
				del self.request.session['alt_files']
			if 'ref_file' in self.request.session.keys():
				del self.request.session['ref_file']
		context['selection'] = {}
		for selection_box, include in self.selection_boxes.items():
			if include:
				context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

		# get attributes of this class and add them to the context
		context['form_code'] = str(self.form_code)
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		# print(context)
		return context



#Class rendering selection box for sequence segments
class SuperpositionWorkflowSelection(AbsSegmentSelection):

	#Left panel
	step = 2
	number_of_steps = 3

	docs = 'structures.html#structure-superposition'

	#Mid section
	#mid_section = 'segment_selection.html'

	#Right panel
	segment_list = True
	buttons = {
		'continue': {
			'label': 'Superpose proteins',
			'url': '/structure/superposition_workflow_results',
			'color': 'success',
		},
	}
	# OrderedDict to preserve the order of the boxes
	selection_boxes = OrderedDict([('reference', False),
		('targets', False),
		('segments', True),])


	def post (self, request, *args, **kwargs):

		# create full selection and import simple selection (if it exists)
		simple_selection = request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		if 'ref_file' in request.FILES:
			request.session['ref_file'] = request.FILES['ref_file']
		if 'alt_files' in request.FILES:
			request.session['alt_files'] = request.FILES.getlist('alt_files')

		context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)
		context['selection'] = {}
		for selection_box, include in self.selection_boxes.items():
			if include:
				context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]


		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		return render(request, self.template_name, context)

	def get_context_data(self, **kwargs):

		self.buttons = {
			'continue': {
				'label': 'Download substructures',
				'url': '/structure/superposition_workflow_results/custom',
				'color': 'success',
			},
		}
		context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)

		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		context['selection'] = {}
		context['selection']['site_residue_groups'] = selection.site_residue_groups
		context['selection']['active_site_residue_group'] = selection.active_site_residue_group
		for selection_box, include in self.selection_boxes.items():
			if include:
				context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

		# get attributes of this class and add them to the context
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		return context


#Class rendering results from superposition workflow
class SuperpositionWorkflowResults(TemplateView):
	"""
	Select download mode for the superposed structures. Full structures, superposed fragments only, select substructure.
	"""

	template_name = 'common_structural_tools.html'

	#Left panel
	step = 3
	number_of_steps = 3
	title = "SELECT SUBSTRUCTURE"
	description = 'Download the desired substructures.'

	#Mid section
	mid_section = 'superposition_results.html'
	#Buttons - none


	def get_context_data (self, **kwargs):

		context = super(SuperpositionWorkflowResults, self).get_context_data(**kwargs)

		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		if 'ref_file' in self.request.session.keys():
			ref_file = StringIO(self.request.session['ref_file'].file.read().decode('UTF-8'))
		elif selection.reference != []:
			ref_file = StringIO(selection.reference[0].item.get_cleaned_pdb())
		if 'alt_files' in self.request.session.keys():
			alt_files = [StringIO(alt_file.file.read().decode('UTF-8')) for alt_file in self.request.session['alt_files']]
		elif selection.targets != []:
			alt_files = [StringIO(x.item.get_cleaned_pdb()) for x in selection.targets if x.type in ['structure', 'structure_model', 'structure_model_Inactive', 'structure_model_Intermediate', 'structure_model_Active']]

		superposition = ProteinSuperpose(deepcopy(ref_file),alt_files, selection)
		out_structs = superposition.run()
		if 'alt_files' in self.request.session.keys():
			alt_file_names = [x.name for x in self.request.session['alt_files']]
		else:
			alt_file_names = []
			for x in selection.targets:
				if x.type=='structure':
					alt_file_names.append('{}_{}.pdb'.format(x.item.protein_conformation.protein.entry_name, x.item.pdb_code.index))
				elif x.type=='structure_model' or x.type=='structure_model_Inactive' or x.type=='structure_model_Intermediate' or x.type=='structure_model_Active':
					alt_file_names.append('Class{}_{}_{}_{}_GPCRdb.pdb'.format(class_dict[x.item.protein.family.slug[:3]], x.item.protein.entry_name, x.item.state.name, x.item.main_template.pdb_code.index))
		if len(out_structs) == 0:
			self.success = False
		elif len(out_structs) >= 1:
			io = PDBIO()
			self.request.session['alt_structs'] = {}
			for alt_struct, alt_file_name in zip(out_structs, alt_file_names):
				tmp = StringIO()
				io.set_structure(alt_struct)
				io.save(tmp)
				self.request.session['alt_structs'][alt_file_name] = tmp

			self.success = True

		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		return context

class SuperpositionWorkflowDownload(View):
	"""
	Serve the (sub)structures depending on user's choice.
	"""

	def get(self, request, *args, **kwargs):

		if self.kwargs['substructure'] == 'select':
			return HttpResponseRedirect('/structure/superposition_workflow_selection')

		io = PDBIO()
		out_stream = BytesIO()
		zipf = zipfile.ZipFile(out_stream, 'w')
		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)
		self.alt_substructure_mapping = {}
		#reference
		if 'ref_file' in request.session.keys():
			self.request.session['ref_file'].file.seek(0)
			ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(self.request.session['ref_file'].file.read().decode('UTF-8')))[0]
			gn_assigner = GenericNumbering(structure=ref_struct)
			gn_assigner.assign_generic_numbers()
			self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
			ref_name = self.request.session['ref_file'].name
		elif selection.reference != []:
			ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(selection.reference[0].item.get_cleaned_pdb()))[0]
			gn_assigner = GenericNumbering(structure=ref_struct)
			gn_assigner.assign_generic_numbers()
			self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
			if selection.reference[0].type=='structure':
				ref_name = '{}_{}_ref.pdb'.format(selection.reference[0].item.protein_conformation.protein.entry_name, selection.reference[0].item.pdb_code.index)
			elif selection.reference[0].type=='structure_model' or selection.reference[0].type=='structure_model_Inactive' or selection.reference[0].type=='structure_model_Intermediate' or selection.reference[0].type=='structure_model_Active':
				ref_name = 'Class{}_{}_{}_{}_GPCRdb_ref.pdb'.format(class_dict[selection.reference[0].item.protein.family.slug[:3]], selection.reference[0].item.protein.entry_name,
																	selection.reference[0].item.state.name, selection.reference[0].item.main_template.pdb_code.index)

		alt_structs = {}
		for alt_id, st in self.request.session['alt_structs'].items():
			st.seek(0)
			alt_structs[alt_id] = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(alt_id, st)[0]
			gn_assigner = GenericNumbering(structure=alt_structs[alt_id])
			gn_assigner.assign_generic_numbers()
			self.alt_substructure_mapping[alt_id] = gn_assigner.get_substructure_mapping_dict()

		if self.kwargs['substructure'] == 'full':

			io.set_structure(ref_struct)
			tmp = StringIO()
			io.save(tmp)
			zipf.writestr(ref_name, tmp.getvalue())

			for alt_name in self.request.session['alt_structs']:
				tmp = StringIO()
				io.set_structure(alt_structs[alt_name])
				io.save(tmp)
				zipf.writestr(alt_name, tmp.getvalue())

		elif self.kwargs['substructure'] == 'substr':

			consensus_gn_set = CASelector(SelectionParser(selection), ref_struct, alt_structs.values()).get_consensus_gn_set()
			io.set_structure(ref_struct)
			tmp = StringIO()
			io.save(tmp, GenericNumbersSelector(consensus_gn_set))
			zipf.writestr(ref_name, tmp.getvalue())
			for alt_name in self.request.session['alt_structs']:
				tmp = StringIO()
				io.set_structure(alt_structs[alt_name])
				io.save(tmp, GenericNumbersSelector(consensus_gn_set))
				zipf.writestr(alt_name, tmp.getvalue())

		elif self.kwargs['substructure'] == 'custom':

			io.set_structure(ref_struct)
			tmp = StringIO()
			io.save(tmp, SubstructureSelector(self.ref_substructure_mapping, parsed_selection=SelectionParser(selection)))

			zipf.writestr(ref_name, tmp.getvalue())
			for alt_name in self.request.session['alt_structs']:
				tmp = StringIO()
				io.set_structure(alt_structs[alt_name])
				io.save(tmp, SubstructureSelector(self.alt_substructure_mapping[alt_name], parsed_selection=SelectionParser(selection)))
				zipf.writestr(alt_name, tmp.getvalue())

		zipf.close()
		if len(out_stream.getvalue()) > 0:
			response = HttpResponse(content_type="application/zip")
			response['Content-Disposition'] = 'attachment; filename="Superposed_structures.zip"'
			response.write(out_stream.getvalue())

		if 'ref_file' in request.FILES:
			request.session['ref_file'] = request.FILES['ref_file']
		if 'alt_files' in request.FILES:
			request.session['alt_files'] = request.FILES.getlist('alt_files')

		return response

# DEPRECATED FUNCTION
class FragmentSuperpositionIndex(TemplateView):

	template_name = 'common_structural_tools.html'

	#Left panel
	step = 1
	number_of_steps = 1

	documentation_url = settings.DOCUMENTATION_URL
	docs = 'sites.html#pharmacophore-generation'

	title = "SUPERPOSE FRAGMENTS OF CRYSTAL STRUCTURES"
	description = """
	The tool implements a fragment-based pharmacophore method, as published in <a href='http://www.ncbi.nlm.nih.gov/pubmed/25286328'>Fidom K, et al (2015)</a>. Interacting ligand moiety - residue pairs extracted from selected crystal structures of GPCRs are superposed onto the input pdb file based on gpcrdb generic residue numbers. Resulting aligned ligand fragments can be used for placement of pharmacophore features.

	Upload a pdb file you want to superpose the interacting moiety - residue pairs.

	Once you have selected all your targets, click the green button.
		"""

	#Input file form data
	header = "Select a file to upload:"
	#Can't control the class properly - staying with the dirty explicit html code
	form_id='fragments'
	form_code = """
	Pdb file:<input id="id_pdb_file" name="pdb_file" type="file" /></br>
	Similarity:</br>
	<input id="similarity" name="similarity" type="radio" value="identical" /> Use fragments with identical residues</br>
	<input checked="checked" id="similarity" name="similarity" type="radio" value="similar" /> Use fragments with residues of similar properties</br>

	Fragments:</br>
	<input checked="checked" id="representative" name="representative" type="radio" value="closest" /> Use fragments from the evolutionary closest crystal structure</br>
	<input id="representative" name="representative" type="radio" value="any" /> Use all available fragments</br></br>
	State:<select id="id_state" name="state">
	<option value="active">Antagonist-bound structures</option>
	<option value="inactive">Agonist-bound structures</option>
	</select>
	"""
	url = '/structure/fragment_superposition_results'
	mid_section = "upload_file_form.html"
	form_height = 350

	#Buttons
	buttons = {
		'continue' : {
			'label' : 'Retrieve fragments',
			'color' : 'success',
			},
		}


	def get_context_data (self, **kwargs):

		context = super(FragmentSuperpositionIndex, self).get_context_data(**kwargs)
		# get attributes of this class and add them to the context
		context['form_code'] = str(self.form_code)
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		return context


# DEPRECATED FUNCTION
class FragmentSuperpositionResults(TemplateView):

	template_name = "common_structural_tools.html"

	#Left panel - blank
	#Mid section
	mid_section = 'fragment_superposition_results.html'
	#Buttons - none

	def post (self, request, *args, **kwargs):

		frag_sp = FragmentSuperpose(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8', 'ignore')),request.FILES['pdb_file'].name)
		superposed_fragments = []
		superposed_fragments_repr = []
		if request.POST['similarity'] == 'identical':
			if request.POST['representative'] == 'any':
				superposed_fragments = frag_sp.superpose_fragments()
			else:
				superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, state=request.POST['state'])
				superposed_fragments = frag_sp.superpose_fragments()
		else:
			if request.POST['representative'] == 'any':
				superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
			else:
				superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, use_similar=True, state=request.POST['state'])
				superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
		if superposed_fragments == []  and superposed_fragments_repr == []:
			self.message = "No fragments were aligned."
		else:
			io = PDBIO()
			out_stream = BytesIO()
			zipf = zipfile.ZipFile(out_stream, 'a', zipfile.ZIP_DEFLATED)
			for fragment, pdb_data in superposed_fragments:
				io.set_structure(pdb_data)
				tmp = StringIO()
				io.save(tmp)
				if request.POST['representative'] == 'any':
					zipf.writestr(fragment.generate_filename(), tmp.getvalue())
				else:
					zipf.writestr("all_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
			if superposed_fragments_repr != []:
				for fragment, pdb_data in superposed_fragments_repr:
					io.set_structure(pdb_data)
					tmp = StringIO()
					io.save(tmp)
					zipf.writestr("representative_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
			zipf.close()
			if len(out_stream.getvalue()) > 0:
				request.session['outfile'] = { 'interacting_moiety_residue_fragments.zip' : out_stream, }
				self.outfile = 'interacting_moiety_residue_fragments.zip'
				self.success = True
				self.zip = 'zip'
				self.message = '{:n} fragments were superposed.'.format(len(superposed_fragments))

		context = super(FragmentSuperpositionResults, self).get_context_data(**kwargs)
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		return render(request, self.template_name, context)


#==============================================================================
class TemplateTargetSelection(AbsReferenceSelection):
	"""
	Starting point for template selection workflow. Target selection.
	"""

	type_of_selection = 'reference'
	# Left panel
	description = 'Select a reference target by searching or browsing in the middle column.' \
		+ '\n\nThe selected reference target will appear in the right column.' \
		+ '\n\nOnce you have selected your reference target, either proceed with all TMs alignment ("Find template"' \
		+ 'button) or specify the sequence segments manualy ("Advanced segment selection" button).'
	step = 1
	number_of_steps = 2
	redirect_on_select = False

	docs = 'structures.html#template-selection'

	# Mid section

	# Right panel
	buttons = OrderedDict()
	buttons['continue'] = {
		'label': 'Find template',
		'url': '/structure/template_browser',
		'color': 'success',
	}
	buttons['segments'] = {
		'label' : 'Advanced segment selection',
		'url' : '/structure/template_segment_selection',
		'color' : 'info',
	}

	selection_boxes = OrderedDict([('reference', True),
		('targets', False),
		('segments', False),])



#==============================================================================
class TemplateSegmentSelection(AbsSegmentSelection):
	"""
	Advanced selection of sequence segments for template search.
	"""
   #Left panel
	step = 2
	number_of_steps = 2

	docs = 'structures.html#template-selection'

	#Mid section
	#mid_section = 'segment_selection.html'

	#Right panel
	segment_list = True
	buttons = {
		'continue': {
			'label': 'Find template',
			'url': '/structure/template_browser',
			'color': 'success',
		},
	}
	# OrderedDict to preserve the order of the boxes
	selection_boxes = OrderedDict([('reference', True),
		('targets', False),
		('segments', True),])



#==============================================================================
class TemplateBrowser(TemplateView):
	"""
	Fetching Structure data and ordering by similarity
	"""

	template_name = "template_browser.html"

	def get_context_data (self, **kwargs):

		context = super(TemplateBrowser, self).get_context_data(**kwargs)

		# get simple selection from session
		simple_selection = self.request.session.get('selection', False)

		# make an alignment
		a = Alignment()
		a.ignore_alternative_residue_numbering_schemes = True

		# load the selected reference into the alignment
		a.load_reference_protein_from_selection(simple_selection)

		# fetch
		qs = Structure.objects.all().select_related(
			"pdb_code__web_resource",
			"protein_conformation__protein__species",
			"protein_conformation__protein__source",
			"protein_conformation__protein__family__parent__parent__parent",
			"publication__web_link__web_resource").prefetch_related(
			"stabilizing_agents",
			"protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
			Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
			annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role')))

		# Dirty but fast
		qsd = {}
		for st in qs:
			qsd[st.protein_conformation.protein.id] = st

		# add proteins to the alignment
		a.load_proteins([x.protein_conformation.protein for x in qs])

		if simple_selection.segments != []:
			a.load_segments_from_selection(simple_selection)
		else:
			a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))

		a.build_alignment()
		a.calculate_similarity()

		context['structures'] = []
		for prot in a.proteins[1:]:
			try:
				context['structures'].append([prot.similarity, prot.identity, qsd[prot.protein.id]])
				del qsd[prot.protein.id]
			except KeyError:
				pass
		return context

#==============================================================================
class PDBClean(TemplateView):
	"""
	Extraction, packing and serving out the pdb records selected via structure/template browser.
	"""
	template_name = "pdb_download.html"

	def get(self, request, *args, **kwargs):
		context = super(PDBClean, self).get_context_data(**kwargs)
		if request.path.endswith('pdb_download_custom'):
			context['trigger_download'] = True

		self.success = False
		self.posted = False
		context['pref'] = True

		# get simple selection from session
		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		if len(selection.targets)>100:
			return HttpResponse("Cannot process more than 100 structures", status=400)
		elif selection.targets != []:
			self.success = True
			context['targets'] = selection.targets

		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		return render(request, self.template_name, context)

	def post(self, request, *args, **kwargs):
		context = super(PDBClean, self).get_context_data(**kwargs)

		self.posted = True
		pref, context['pref'] = True, True
		water, context['water'] = False, False
		hets, context['hets'] = False, False

		if 'pref_chain' not in request.POST.keys():
			pref, context['pref'] = False, False
		if 'water' in request.POST.keys():
			water, context['water'] = True, True
		if 'hets' in request.POST.keys():
			hets, context['hets'] = True, True

		# get simple selection from session
		simple_selection = request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)
		out_stream = BytesIO()
		io = PDBIO()
		zipf = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
		structs = []
		if selection.targets != []:
			if selection.targets != [] and selection.targets[0].type == 'structure':
				request.session['substructure_mapping'] = OrderedDict()
				for selected_struct in [x for x in selection.targets if x.type == 'structure']:
					struct_name = '{}_{}.pdb'.format(selected_struct.item.protein_conformation.protein.parent.entry_name, selected_struct.item.pdb_code.index)
					if hets:
						lig_names = [x.pdb_reference for x in StructureLigandInteraction.objects.filter(structure=selected_struct.item, annotated=True)]
					else:
						lig_names = None
					gn_assigner = GenericNumberingFromDB(selected_struct.item, PDBParser(QUIET=True).get_structure(struct_name, StringIO(selected_struct.item.get_cleaned_pdb(pref, water, lig_names)))[0])
					tmp = StringIO()
					io.set_structure(gn_assigner.assign_generic_numbers())
					request.session['substructure_mapping'][struct_name] = gn_assigner.get_substructure_mapping_dict()
					io.save(tmp)
					zipf.writestr(struct_name, tmp.getvalue())
					del gn_assigner, tmp
					structs.append(selected_struct.item.id)
				for struct in selection.targets:
					selection.remove('targets', 'structure', struct.item.id)
				for struct_id in structs:
					s = Structure.objects.get(id=struct_id)
					selection.add('targets', 'structure', SelectionItem('structure', s))

			elif selection.targets != [] and selection.targets[0].type in ['structure_model', 'structure_model_Inactive', 'structure_model_Intermediate', 'structure_model_Active']:
				for hommod in [x for x in selection.targets if x.type in ['structure_model', 'structure_model_Inactive', 'structure_model_Intermediate', 'structure_model_Active']]:
					mod_name = 'Class{}_{}_{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.item.protein.family.slug[:3]], hommod.item.protein.entry_name,
																				  hommod.item.state.name, hommod.item.main_template.pdb_code.index, hommod.item.version)
					tmp = StringIO(hommod.item.pdb)
					request.session['substructure_mapping'] = 'full'
					zipf.writestr(mod_name, tmp.getvalue())
					del tmp

					# stat file
					# rotamers = StructureModelStatsRotamer.objects.filter(homology_model=hommod.item).prefetch_related('residue','backbone_template','rotamer_template').order_by('residue__sequence_number')
					# stats_data = 'Segment,Sequence_number,Generic_number,Backbone_template,Rotamer_template\n'
					# for r in rotamers:
					#     try:
					#         gn = r.residue.generic_number.label
					#     except:
					#         gn = '-'
					#     if r.backbone_template:
					#         bt = r.backbone_template.pdb_code.index
					#     else:
					#         bt = '-'
					#     if r.rotamer_template:
					#         rt = r.rotamer_template.pdb_code.index
					#     else:
					#         rt = '-'
					#     stats_data+='{},{},{},{},{}\n'.format(r.residue.protein_segment.slug, r.residue.sequence_number, gn, bt, rt)
					# stats_name = mod_name[:-3]+'templates.csv'
					# zipf.writestr(stats_name, stats_data)
					# del stats_data

				for mod in selection.targets:
					selection.remove('targets', 'structure_model', mod.item.id)

			# export simple selection that can be serialized
			simple_selection = selection.exporter()

			request.session['selection'] = simple_selection
			request.session['cleaned_structures'] = out_stream
		self.success = True
		context['prepared_structures'] = True
		context['targets'] = selection.targets
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]

		if selection.targets != [] and selection.targets[0].type == 'structure_model':
			return zipf
		else:
			return render(request, self.template_name, context)


class PDBSegmentSelection(AbsSegmentSelection):

	#Left panel
	step = 2
	number_of_steps = 2

	#Mid section
	#mid_section = 'segment_selection.html'

	#Right panel
	segment_list = True

	# OrderedDict to preserve the order of the boxes
	selection_boxes = OrderedDict([('reference', False),
		('targets', False),
		('segments', True),])

	def get_context_data(self, **kwargs):

		self.buttons = {
			'continue': {
				'label': 'Download substructures',
				'url': '/structure/pdb_download_custom',
				'color': 'success',
			},
		}
		context = super(PDBSegmentSelection, self).get_context_data(**kwargs)

		simple_selection = self.request.session.get('selection', False)
		selection = Selection()
		if simple_selection:
			selection.importer(simple_selection)

		context['selection'] = {}
		# context['selection']['site_residue_groups'] = selection.site_residue_groups
		# context['selection']['active_site_residue_group'] = selection.active_site_residue_group

		context['ignore_residue_selection'] = True
		for selection_box, include in self.selection_boxes.items():
			if include:
				context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

		# get attributes of this class and add them to the context
		attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
		for a in attributes:
			if not(a[0].startswith('__') and a[0].endswith('__')):
				context[a[0]] = a[1]
		return context

#==============================================================================
class PDBDownload(View):
	"""
	Serve the PDB (sub)structures depending on user's choice.
	"""

	def get(self, request, hommods=False, *args, **kwargs):
		if self.kwargs['substructure'] == 'select':
			return HttpResponseRedirect('/structure/pdb_segment_selection')

		if self.kwargs['substructure'] == 'full':
			out_stream = request.session['cleaned_structures']

		elif self.kwargs['substructure'] == 'custom':
			simple_selection = request.session.get('selection', False)
			selection = Selection()
			if simple_selection:
				selection.importer(simple_selection)
			io = PDBIO()
			zipf_in = zipfile.ZipFile(request.session['cleaned_structures'], 'r')
			out_stream = BytesIO()
			zipf_out = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
			for name in zipf_in.namelist():
				tmp = StringIO()
				io.set_structure(PDBParser(QUIET=True).get_structure(name, StringIO(zipf_in.read(name).decode('utf-8')))[0])
				io.save(tmp, SubstructureSelector(request.session['substructure_mapping'][name], parsed_selection=SelectionParser(selection)))
				zipf_out.writestr(name, tmp.getvalue())

			zipf_in.close()
			zipf_out.close()

		if len(out_stream.getvalue()) > 0:
			response = HttpResponse(content_type="application/zip")
			if hommods == False:
				response['Content-Disposition'] = 'attachment; filename="pdb_structures.zip"'
			else:
				response['Content-Disposition'] = 'attachment; filename="GPCRDB_homology_models.zip"'
			response.write(out_stream.getvalue())

		return response

#==============================================================================
def ConvertStructuresToProteins(request):
	"For alignment from structure browser"

	simple_selection = request.session.get('selection', False)
	selection = Selection()
	if simple_selection:
		selection.importer(simple_selection)
	if selection.targets != []:
		for struct in selection.targets:
			prot = struct.item.protein_conformation.protein.parent
			selection.remove('targets', 'structure', struct.item.id)
			selection.add('targets', 'protein', SelectionItem('protein', prot))
		if selection.reference != []:
			selection.add('targets', 'protein', selection.reference[0])
	# export simple selection that can be serialized
	simple_selection = selection.exporter()

	# add simple selection to session
	request.session['selection'] = simple_selection

	return HttpResponseRedirect('/alignment/segmentselection')


def ConvertStructureModelsToProteins(request):
	"For alignment from homology model browser"
	simple_selection = request.session.get('selection', False)
	selection = Selection()
	if simple_selection:
		selection.importer(simple_selection)
	if selection.targets != []:
		for struct_mod in selection.targets:
			if hasattr(struct_mod.item, 'protein'):
				if not struct_mod.item.protein.accession:
					prot = struct_mod.item.protein.parent
				else:
					prot = struct_mod.item.protein
				selection.remove('targets', 'structure_model', struct_mod.item.id)
				selection.add('targets', 'protein', SelectionItem('protein', prot))
			elif hasattr(struct_mod.item, 'receptor_protein'):
				if not struct_mod.item.receptor_protein.accession:
					prot = struct_mod.item.receptor_protein.parent
				else:
					prot = struct_mod.item.receptor_protein
				selection.remove('targets', 'structure_complex_receptor', struct_mod.item.id)
				selection.add('targets', 'protein', SelectionItem('protein', prot))
			elif hasattr(struct_mod.item, 'pdb_code'):
				prot = struct_mod.item.protein_conformation.protein.parent
				selection.remove('targets', 'structure', struct_mod.item.id)
				selection.add('targets', 'protein', SelectionItem('protein', prot))

		if selection.reference != []:
			selection.add('targets', 'protein', selection.reference[0])
	# export simple selection that can be serialized
	simple_selection = selection.exporter()

	# add simple selection to session
	request.session['selection'] = simple_selection

	return HttpResponseRedirect('/alignment/segmentselection')

def ConvertStructureComplexSignprotToProteins(request):
	"For alignment from complex model browser, specifically for signprots"
	simple_selection = request.session.get('selection', False)
	selection = Selection()
	if simple_selection:
		selection.importer(simple_selection)
	if selection.targets != []:
		for struct_mod in selection.targets:
			if hasattr(struct_mod.item, 'sign_protein'):
				prot = struct_mod.item.sign_protein
				selection.remove('targets', 'structure_complex_signprot', struct_mod.item.id)
				selection.add('targets', 'protein', SelectionItem('protein', prot))
			else:
				prot = SignprotComplex.objects.get(structure=struct_mod.item).protein
				selection.remove('targets', 'structure', struct_mod.item.id)
				selection.add('targets', 'protein', SelectionItem('protein', prot))
		if selection.reference != []:
			selection.add('targets', 'protein', selection.reference[0])
	# export simple selection that can be serialized
	simple_selection = selection.exporter()

	# add simple selection to session
	request.session['selection'] = simple_selection

	return HttpResponseRedirect('/alignment/segmentselectiongprot')


def HommodDownload(request):
	"Download selected homology models in zip file"
	pks = request.GET['ids'].split(',')

	hommodels = [StructureModel.objects.get(pk=pk) for pk in pks]

	zip_io = BytesIO()
	with zipfile.ZipFile(zip_io, mode='w', compression=zipfile.ZIP_DEFLATED) as backup_zip:
		for hommod in hommodels:
			io = StringIO(hommod.pdb_data.pdb)
			stats_text = StringIO(hommod.stats_text.stats_text)
			if not hommod.protein.accession:
				mod_name = 'Class{}_{}_{}_refined_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.parent.entry_name,
																   hommod.main_template.pdb_code.index, hommod.state.name, hommod.version)
				stat_name = 'Class{}_{}_{}_refined_{}_{}_GPCRDB.templates.csv'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.parent.entry_name,
																   hommod.main_template.pdb_code.index, hommod.state.name, hommod.version)
			else:
				mod_name = 'Class{}_{}_{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name,
																		  hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
				stat_name = 'Class{}_{}_{}_{}_{}_GPCRDB.templates.csv'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name,
																		  hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
			backup_zip.writestr(mod_name, io.getvalue())
			backup_zip.writestr(stat_name, stats_text.getvalue())

	response = HttpResponse(zip_io.getvalue(), content_type='application/x-zip-compressed')
	response['Content-Disposition'] = 'attachment; filename=%s' % 'GPCRDB_homology_models' + ".zip"
	response['Content-Length'] = zip_io.tell()
	return response

def ComplexmodDownload(request):
	"Download selected complex homology models in zip file"
	pks = request.GET['ids'].split(',')

	hommodels = StructureComplexModel.objects.filter(pk__in=pks).prefetch_related('receptor_protein__family','main_template__pdb_code').all()
	zip_io = BytesIO()
	with zipfile.ZipFile(zip_io, mode='w', compression=zipfile.ZIP_DEFLATED) as backup_zip:
		for hommod in hommodels:
			io = StringIO(hommod.pdb_data.pdb)
			stats_text = StringIO(hommod.stats_text.stats_text)
			if not hommod.receptor_protein.accession:
				mod_name = 'Class{}_{}-{}_{}_refined_{}_GPCRDB.pdb'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.parent.entry_name,
																   hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
				stat_name = 'Class{}_{}-{}_{}_refined_{}_GPCRDB.templates.csv'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.parent.entry_name,
																   hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
			else:
				mod_name = 'Class{}_{}-{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.entry_name,
																   hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
				stat_name = 'Class{}_{}-{}_{}_{}_GPCRDB.templates.csv'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.entry_name,
																   hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
			backup_zip.writestr(mod_name, io.getvalue())
			backup_zip.writestr(stat_name, stats_text.getvalue())
	response = HttpResponse(zip_io.getvalue(), content_type='application/x-zip-compressed')
	response['Content-Disposition'] = 'attachment; filename=%s' % 'GPCRDB_complex_homology_models' + ".zip"
	response['Content-Length'] = zip_io.tell()
	return response

def SingleModelDownload(request, modelname, fullness, state=None, csv=False):
	"Download single homology model"
	zip_io = BytesIO()
	if state:
		hommod = StructureModel.objects.get(protein__entry_name=modelname.lower(), state__slug=state)
	else:
		hommod = StructureModel.objects.get(protein__entry_name=modelname.lower())
	if not hommod.protein.accession:
		version = hommod.pdb_data.pdb.split('\n')[0][-10:]
		mod_name = 'Class{}_{}_{}_refined_{}_{}_GPCRdb.pdb'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.parent.entry_name,
																 hommod.main_template.pdb_code.index, hommod.state.name, version)
		stat_name = 'Class{}_{}_{}_refined_{}_{}_GPCRdb.templates.csv'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.parent.entry_name,
																 hommod.main_template.pdb_code.index, hommod.state.name, version)
	else:
		mod_name = 'Class{}_{}_{}_{}_{}_GPCRdb.pdb'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name,
																	   hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
		stat_name = 'Class{}_{}_{}_{}_{}_GPCRdb.templates.csv'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name,
																	   hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
	pdb_lines = hommod.pdb_data.pdb
	stats_lines = hommod.stats_text.stats_text
	if fullness=='noloops':
		filtered_lines = ''
		remark_line_count = 0
		for line in pdb_lines.split('\n'):
			if line.startswith('REMARK'):
				filtered_lines+=line+'\n'
				remark_line_count+=1
		filtered_lines+='REMARK    {} ALL LOOPS REMOVED AFTER MODELING\n'.format(remark_line_count+1)
		if not hommod.protein.accession:
			helix_resis = Residue.objects.filter(protein_conformation__protein=hommod.protein.parent, protein_segment__category='helix').values_list('sequence_number', flat=True)
		else:
			helix_resis = Residue.objects.filter(protein_conformation__protein=hommod.protein, protein_segment__category='helix').values_list('sequence_number', flat=True)
		p = PDBParser(get_header=True)
		pdb = p.get_structure('pdb', StringIO(pdb_lines))[0]
		to_remove = []
		for chain in pdb:
			for res in chain:
				if res.id[1] not in helix_resis and res.id[0]==' ':
					to_remove.append(res.id)
		for i in to_remove:
			chain.detach_child(i)
		io = StringIO()
		pdbio = PDBIO()
		pdbio.set_structure(pdb)
		pdbio.save(io)
		io = StringIO(filtered_lines+io.getvalue())
		filtered_stats_lines = ''
		for i in stats_lines.split('\n'):
			if i.startswith('Segment'):
				filtered_stats_lines+=i+'\n'
				continue
			if len(i)<2:
				continue
			resnum = int(i.split(',')[1])
			if resnum in helix_resis:
				filtered_stats_lines+=i+'\n'
		stats_text = StringIO(filtered_stats_lines)
		mod_name = mod_name.split('.')[0]+'_WL'+'.pdb'
		stat_name = stat_name.split('.')[0]+'_WL'+'.templates.csv'
	else:
		io = StringIO(pdb_lines)
		stats_text = StringIO(stats_lines)
	with zipfile.ZipFile(zip_io, mode='w', compression=zipfile.ZIP_DEFLATED) as backup_zip:
		backup_zip.writestr(mod_name, io.getvalue())
		backup_zip.writestr(stat_name, stats_text.getvalue())
	response = HttpResponse(zip_io.getvalue(), content_type='application/x-zip-compressed')
	response['Content-Disposition'] = 'attachment; filename=%s' % mod_name.split('.')[0] + ".zip"
	response['Content-Length'] = zip_io.tell()

	return response

def SingleComplexModelDownload(request, modelname, signprot, csv=False):
	"Download single homology model"

	zip_io = BytesIO()
	if signprot=='complex':
		hommod = StructureComplexModel.objects.get(receptor_protein__entry_name=modelname.lower())
		mod_name = 'Class{}_{}-{}_{}_refined_{}_GPCRdb.pdb'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.parent.entry_name,
															hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
		stat_name = 'Class{}_{}-{}_{}_refined_{}_GPCRdb.templates.csv'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.parent.entry_name,
															hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
	else:
		hommod = StructureComplexModel.objects.get(receptor_protein__entry_name=modelname, sign_protein__entry_name=signprot)
		mod_name = 'Class{}_{}-{}_{}_{}_GPCRdb.pdb'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.entry_name,
															hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
		stat_name = 'Class{}_{}-{}_{}_{}_GPCRdb.templates.csv'.format(class_dict[hommod.receptor_protein.family.slug[:3]], hommod.receptor_protein.entry_name,
															hommod.sign_protein.entry_name, hommod.main_template.pdb_code.index, hommod.version)
	io = StringIO(hommod.pdb_data.pdb)
	stats_text = StringIO(hommod.stats_text.stats_text)
	with zipfile.ZipFile(zip_io, mode='w', compression=zipfile.ZIP_DEFLATED) as backup_zip:
		backup_zip.writestr(mod_name, io.getvalue())
		backup_zip.writestr(stat_name, stats_text.getvalue())
	response = HttpResponse(zip_io.getvalue(), content_type='application/x-zip-compressed')
	response['Content-Disposition'] = 'attachment; filename=%s' % mod_name.split('.')[0] + ".zip"
	response['Content-Length'] = zip_io.tell()

	return response

def ServePdbOutfile (request, outfile, replacement_tag):

	root, ext = os.path.splitext(outfile)
	out_stream = request.session['outfile'][outfile]
	response = HttpResponse(content_type="chemical/x-pdb")
	response['Content-Disposition'] = 'attachment; filename="{}_{}.pdb"'.format(root, replacement_tag)
	response.write(out_stream.getvalue())

	return response


def ServeZipOutfile (request, outfile):

	out_stream = request.session['outfile'][outfile]
	response = HttpResponse(content_type="application/zip")
	response['Content-Disposition'] = 'attachment; filename="{}"'.format(outfile)
	response.write(out_stream.getvalue())

	return response

def RenderTrees(request):
	number = request.GET['number']
	tree = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_tree.xml').read()
	legend = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_legend.svg').read()
	context = {'tree':tree, 'leg':legend, 'num':number}
	return render(request, 'phylogenetic_trees.html', context)

def webform(request):
	form = construct_form()
	context = {'form':form}
	return render(request, 'web_form.html',context)

def webform_two(request, slug=None):
	context = {}
	if slug:
		c = Construct.objects.filter(name=slug).get()
		# print(c.json)
		# test = ast.literal_eval(c.json)
		# print(test)
		json_data = json.loads(c.json)
		if 'raw_data' not in json_data:
			json_data = convert_ordered_to_disordered_annotation(json_data)
		else:
			if 'csrfmiddlewaretoken' in json_data['raw_data']:
				del json_data['raw_data']['csrfmiddlewaretoken'] #remove to prevent errors

		context = {'edit':json.dumps(json_data)}
	return render(request, 'web_form_2.html',context)

def webformdata(request) :

	data = request.POST
	raw_data = deepcopy(data)
	purge_keys = ('Please Select','aamod_position','wt_aa','mut_aa','insert_pos_type','protein_type','deletion','csrfmiddlewaretoken')
	data = dict((k, v) for k, v in data.items() if v!='' and v!='Please Select') #remove empty
	deletions = []
	mutations = []
	contact_info= OrderedDict()
	construct_crystal=OrderedDict()
	auxiliary=OrderedDict()
	expression=OrderedDict()
	solubilization = OrderedDict()
	crystallization = OrderedDict()
	modifications = []

	error = 0
	error_msg = []
	for key,value in sorted(data.items()):
		try:
			if key.startswith('delet_start'):
				deletions.append({'start':value, 'end':data[key.replace('start','end')], 'origin':'user', 'type':'range'})
				data.pop(key, None)
				data.pop(key.replace('start','end'), None)
			elif key.startswith('ins_start'):
				deletions.append({'start':value, 'end':data[key.replace('start','end')], 'origin':'insertion'+key.replace('ins_start',''), 'type':'range'})
				data.pop(key, None)
				data.pop(key.replace('start','end'), None)
				data.pop(key.replace('ins_start',''), None)
			elif key.startswith(('deletion_single', 'insert_pos_single')):
				if key.startswith('insert_pos_single'):
					deletions.append({'pos':value, 'origin':'insertion'+key.replace('insert_pos_single',''), 'type':'single'})
					data.pop(key.replace('insert_pos_single',''), None)
				else:
					deletions.append({'pos':value, 'origin':'user', 'type':'single'})
				data.pop(key, None)

			if key.startswith('aa_no'):
				pos_id = key.replace('aa_no','')
				# if pos_id=='':
				# 	mut_id='1'
				# else:
				# 	mut_id=pos_id.replace('_','')

				if 'mut_remark'+pos_id in data:
					remark = data['mut_remark'+pos_id]
				else:
					remark = ''

				mutations.append({'pos':value,'wt':data['wt_aa'+pos_id],'mut':data['mut_aa'+pos_id], 'type':data['mut_type'+pos_id], 'remark':remark})
				data.pop(key, None)
				data.pop('wt_aa'+pos_id, None)
				data.pop('mut_aa'+pos_id, None)
				data.pop('mut_type'+pos_id, None)

			if key.startswith(('date','name_cont', 'pi_name',
				'pi_address','address','url','pi_email' )):
				contact_info[key]=value
				data.pop(key, None)

			if key.startswith(('pdb', 'pdb_name',
				'uniprot','ligand_name', 'ligand_activity', 'ligand_conc', 'ligand_conc_unit','ligand_id','ligand_id_type')):
				construct_crystal[key]=value
				data.pop(key, None)

			if key.startswith('position'):
				pos_id = key.replace('position','')
				if pos_id=='':
					aux_id='1'
				else:
					aux_id=pos_id.replace('_','')

				if 'aux'+aux_id not in auxiliary:
					auxiliary['aux'+aux_id] = {'position':value,'type':data['protein_type'+pos_id],'presence':data['presence'+pos_id]}

					data.pop(key, None)
					data.pop('protein_type'+pos_id, None)
					data.pop('presence'+pos_id, None)

			if key.startswith(('tag', 'fusion_prot', 'signal', 'linker_seq','prot_cleavage', 'other_prot_cleavage' )):
				temp = key.split('_')
				if len(temp)==4:
					pos_id = "_"+temp[3]
					aux_id=pos_id.replace('_','')
				elif len(temp)==3:
					pos_id = "_"+temp[2]
					aux_id=pos_id.replace('_','')
				elif len(temp)==2 and temp[1].isdigit():
					pos_id = "_"+temp[1]
					aux_id=pos_id.replace('_','')
				else:
					pos_id = ''
					aux_id = '1'
				# print(key,aux_id,pos_id)

				if 'aux'+aux_id not in auxiliary:
					auxiliary['aux'+aux_id] = {'position':data['position'+pos_id],'type':data['protein_type'+pos_id],'presence':data['presence'+pos_id]}

					data.pop('position'+pos_id, None)
					data.pop('protein_type'+pos_id, None)
					data.pop('presence'+pos_id, None)

				# if value=='Other':
				#     auxiliary['aux'+aux_id]['other'] = data['other_'+auxiliary['aux'+aux_id]['type']+pos_id]
				#     data.pop('other_'+auxiliary['aux'+aux_id]['type']+pos_id,None)

				auxiliary['aux'+aux_id]['subtype'] = value
				data.pop(key, None)

			if key.startswith(('expr_method', 'host_cell_type',
					'host_cell', 'expr_remark','expr_other','other_host','other_host_cell' )):
				expression[key]=value
				data.pop(key, None)

			if key.startswith(('deterg_type','deterg_concentr','deterg_concentr_unit','solub_additive','additive_concentr','addit_concentr_unit','chem_enz_treatment','sol_remark')):
				solubilization[key]=value
				data.pop(key, None)

			elif key.startswith(('crystal_type','crystal_method','other_method','other_crystal_type',
							   'protein_concentr','protein_conc_unit','temperature','ph_single','ph',
							   'ph_range_one','ph_range_two','crystal_remark','lcp_lipid','lcp_add',
							   'lcp_conc','lcp_conc_unit','detergent','deterg_conc','deterg_conc_unit','lipid','lipid_concentr','lipid_concentr_unit',
							   'other_deterg','other_deterg_type', 'other_lcp_lipid','other_lipid')):
				crystallization[key]=value
				data.pop(key, None)

			if key.startswith('chemical_comp') and not key.startswith('chemical_comp_type'):

				if 'chemical_components' not in crystallization:
					crystallization['chemical_components'] = []

				# print(key)
				if key!='chemical_comp': #not first
					comp_id = key.replace('chemical_comp','')
				else:
					comp_id = ''

				crystallization['chemical_components'].append({'component':value,'type':data['chemical_comp_type'+comp_id],'value':data['concentr'+comp_id],'unit':data['concentr_unit'+comp_id]})
				data.pop(key, None)
				data.pop('concentr'+comp_id, None)
				data.pop('concentr_unit'+comp_id, None)
				data.pop('chemical_comp_type'+comp_id, None)


			if key.startswith('aamod') and not key.startswith('aamod_position') and not key.startswith('aamod_pair') and not key=='aamod_position' and not key=='aamod_single':
				if key!='aamod': #not first
					mod_id = key.replace('aamod','')
				else:
					mod_id = ''

				if data['aamod_position'+mod_id]=='single':
					pos = ['single',data['aamod_single'+mod_id]]
					data.pop('aamod_single'+mod_id, None)
				elif data['aamod_position'+mod_id]=='range':
					pos = ['range',[data['aamod_start'+mod_id],data['aamod_end'+mod_id]]]
					data.pop('aamod_start'+mod_id, None)
					data.pop('aamod_end'+mod_id, None)
				elif data['aamod_position'+mod_id]=='pair':
					pos = ['pair',[data['aamod_pair_one'+mod_id],data['aamod_pair_two'+mod_id]]]
					data.pop('aamod_pair_one'+mod_id, None)
					data.pop('aamod_pair_two'+mod_id, None)

				remark = ''
				if 'mod_remark'+mod_id in data:
					remark = data['mod_remark'+mod_id]
				modifications.append({'type':value,'remark':remark,'position':pos })
				data.pop(key, None)
				data.pop('mod_remark'+mod_id, None)
				data.pop('aamod_position'+mod_id, None)

			if key.startswith(purge_keys):
				data.pop(key, None)
		except BaseException as e:
			error_msg.append(str(e))
			error = 1

	auxiliary = OrderedDict(sorted(auxiliary.items()))

	context = OrderedDict( [('contact_info',contact_info), ('construct_crystal',construct_crystal),
						   ('auxiliary' , auxiliary),  ('deletions',deletions), ('mutations',mutations),
						   ('modifications', modifications), ('expression', expression), ('solubilization',solubilization),
						   ('crystallization',crystallization),  ('unparsed',data),  ('raw_data',raw_data), ('error', error), ('error_msg',error_msg)] )

	add_construct(context)

	if error==0:
		dump_dir = '/protwis/construct_dump'
		# dump_dir = '/web/sites/files/construct_data' #for sites
		if not os.path.exists(dump_dir):
			os.makedirs(dump_dir)
		ts = int(time.time())
		json_data = context
		json.dump(json_data, open(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json", 'w'), indent=4, separators=(',', ': '))

		context['data'] = sorted(data.items())
		#context['data'] = sorted(raw_data.items())

		recipients = ['christian@munk.be']
		emaillist = [elem.strip().split(',') for elem in recipients]
		msg = MIMEMultipart()
		msg['Subject'] = 'GPCRdb: New webform data'
		msg['From'] = 'gpcrdb@gmail.com'
		msg['Reply-to'] = 'gpcrdb@gmail.com'

		msg.preamble = 'Multipart massage.\n'

		part = MIMEText("Hi, please find the attached file")
		msg.attach(part)

		part = MIMEApplication(open(str(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json"),"rb").read())
		part.add_header('Content-Disposition', 'attachment', filename=str(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json"))
		msg.attach(part)


		server = smtplib.SMTP("smtp.gmail.com:587")
		server.ehlo()
		server.starttls()
		server.login("gpcrdb@gmail.com", "gpcrdb2016")

		server.sendmail(msg['From'], emaillist , msg.as_string())

		context['filename'] = str(ts)+"_"+construct_crystal['pdb']

	return render(request, 'web_form_results.html', context)

def webform_download(request,slug):
	dump_dir = '/protwis/construct_dump'
	# dump_dir = '/web/sites/files/construct_data' #for sites
	file = dump_dir+"/"+str(slug)+".json"
	out_stream = open(file,"rb").read()
	response = HttpResponse(content_type="application/json")
	response['Content-Disposition'] = 'attachment; filename="{}"'.format(file)
	response.write(out_stream)
	return response
