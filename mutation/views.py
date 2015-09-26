from django.shortcuts import render
from django.template import loader, Context

from django.http import HttpResponse
from django.conf import settings
from mutation.functions import *
from mutation.models import *

from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from residue.models import Residue,ResidueNumberingScheme
from residue.views import ResidueTablesDisplay
from protein.models import Protein,ProteinSegment


from datetime import datetime
from collections import OrderedDict
import json
#env/bin/python3 -m pip install xlrd

from io import BytesIO
import re
import math
import urllib
import xlsxwriter #sudo pip3 install XlsxWriter

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    docs = '/docs/alignment'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/mutations/segmentselection',
            'color': 'success',
        },
    }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = '/docs/alignment'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show mutants',
            'url': '/mutations/render',
            'color': 'success',
        },
    }

def render_mutations(request, protein = None, family = None, download = None, **response_kwargs):  

    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    
     # local protein list
    proteins = []

    if protein: #if protein static page

        proteins.append(Protein.objects.get(entry_name = protein))
        segments_ids = ProteinSegment.objects.all().values('id')
        original_segments = ProteinSegment.objects.all()

    elif family:

        family_proteins = Protein.objects.filter(family__slug__startswith=family, sequence_type__slug='wt')
        for fp in family_proteins:
            proteins.append(fp)
        segments_ids = ProteinSegment.objects.all().values('id')
        original_segments = ProteinSegment.objects.all()

    else:

        # flatten the selection into individual proteins
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)
            elif target.type == 'family':
                # species filter
                species_list = []
                for species in simple_selection.species:
                    species_list.append(species.item)

                # annotation filter
                protein_source_list = []
                for protein_source in simple_selection.annotation:
                    protein_source_list.append(protein_source.item)
                    
                family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                    species__in=(species_list),
                    source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                for fp in family_proteins:
                    proteins.append(fp)

        original_segments = []
        segments_ids = []
        for segment in simple_selection.segments:
            original_segments.append(segment.item)
            segments_ids.append(segment.item.id)

       #scheme
    used_schemes = {}
    for protein in proteins:
        if protein.residue_numbering_scheme.slug not in used_schemes:
            used_schemes[protein.residue_numbering_scheme.slug] = 0
        used_schemes[protein.residue_numbering_scheme.slug] += 1

    used_scheme = max(used_schemes, key=used_schemes.get)

    #print(segments)

    mutations = MutationExperiment.objects.filter(protein__in=proteins, 
                            residue__protein_segment__in=original_segments).prefetch_related('protein', 
                            'residue__protein_segment','residue__display_generic_number',
                            'residue__generic_number', 'residue','exp_func','exp_qual',
                            'exp_measure', 'exp_type', 'ligand_role', 'ligand','refs','raw',
                            'ligand__properities', 'refs__web_link', 'refs__web_link__web_resource')
    
    mutations_list = [x.residue.display_generic_number.label for x in mutations if x.residue.display_generic_number]
    
    mutations_list = {}
    mutations_display_generic_number = {}
    for mutation in mutations:
        if not mutation.residue.display_generic_number: continue #cant map those without display numbers
        if mutation.residue.display_generic_number.label not in mutations_list: mutations_list[mutation.residue.display_generic_number.label] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        mutations_list[mutation.residue.display_generic_number.label].append([mutation.foldchange,ligand,qual])

        mutations_display_generic_number[mutation.raw.id] = mutation.residue.display_generic_number.label

    print(mutations_display_generic_number)

    # create an alignment object
    a = Alignment()

    a.load_proteins(proteins)
    segments = ProteinSegment.objects.all().filter().prefetch_related()
    a.load_segments(segments) #get all segments to make correct diagrams

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    residue_list = []
    generic_numbers = []
    reference_generic_numbers = {}
    count = 1 #build sequence_number

    ################################################################################
    #FIXME -- getting matching display_generic_numbers kinda randomly.

    for seg in a.consensus: #Grab list of generic_numbers to lookup for their display numbers
        for aa in a.consensus[seg]:
            if "x" in aa:
                generic_numbers.append(aa)

    generic_ids = Residue.objects.filter(generic_number__label__in=generic_numbers,display_generic_number__scheme__slug=used_scheme).values('id').distinct('generic_number__label').order_by('generic_number__label')
    rs = Residue.objects.filter(id__in=generic_ids).prefetch_related('display_generic_number','generic_number')

    for r in rs: #make lookup dic.
        reference_generic_numbers[r.generic_number.label] = r
    ################################################################################
    mutations_pos_list = {}
    for seg in a.consensus:
        for aa,v in a.consensus[seg].items():
            r = Residue()
            r.sequence_number =  count #FIXME is this certain to be correct that the position in consensus is seq position? 
            if "x" in aa:
                r.display_generic_number = reference_generic_numbers[aa].display_generic_number #FIXME
                if r.display_generic_number.label in mutations_list:
                    if r.sequence_number not in mutations_pos_list: 
                        mutations_pos_list[r.sequence_number] = []
                    mutations_pos_list[r.sequence_number].append(mutations_list[r.display_generic_number.label])
                r.segment_slug = seg
                r.family_generic_number = aa
            else:
                r.segment_slug = seg
                r.family_generic_number = aa
            r.amino_acid = v[0]
            r.extra = v[2] #Grab consensus information
            residue_list.append(r)

            count += 1         

    protein_ids = list(set([x.id for x in proteins]))
    HelixBox = DrawHelixBox(residue_list,'Class A',str(protein_ids), nobuttons = 1)
    SnakePlot = DrawSnakePlot(residue_list,'Class A',str(protein_ids), nobuttons = 1)

    context = {}
    numbering_schemes_selection = ['gpcrdba','gpcrdbb','gpcrdbc','gpcrdbf'] #there is a residue_numbering_scheme attribute on the protein model, so it's easy to find out
    numbering_schemes_selection = list(used_schemes.keys())
    numbering_schemes_selection = ['gpcrdb'] + numbering_schemes_selection #always use A for reference
    numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()

    segments = ProteinSegment.objects.filter(pk__in=segments_ids,category='helix')

    if ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME) in numbering_schemes:
        default_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
    else:
        default_scheme = numbering_schemes[0]

    residue_table_list = []
    for mutation in mutations:
        residue_table_list.append(mutation.residue.generic_number)

    # prepare the dictionary
    # each helix has a dictionary of positions
    # default_generic_number or first scheme on the list is the key
    # value is a dictionary of other gn positions and residues from selected proteins 
    data = OrderedDict()
    for segment in segments:
        data[segment.slug] = OrderedDict()
        # residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins).prefetch_related('protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        #     'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
        residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins, 
                                        generic_number__in=residue_table_list).prefetch_related('protein_conformation__protein', 
                                        'protein_conformation__state', 'protein_segment',
                                        'generic_number','display_generic_number','generic_number__scheme', 'alternative_generic_numbers__scheme')
        for scheme in numbering_schemes:
            if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
            elif scheme == default_scheme:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}

        for residue in residues:
            alternatives = residue.alternative_generic_numbers.all()
            pos = residue.generic_number
            for alternative in alternatives:
                if alternative.scheme not in numbering_schemes:
                    continue
                scheme = alternative.scheme
                if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    pos = residue.generic_number
                    if scheme == pos.scheme:
                        data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                    else:
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][scheme.slug] = alternative.label
                        if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                            data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                        data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                else:
                    if scheme.slug not in data[segment.slug][pos.label].keys():
                        data[segment.slug][pos.label][scheme.slug] = alternative.label
                    if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                        data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                    data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)

    # Preparing the dictionary of list of lists. Dealing with tripple nested dictionary in django templates is a nightmare
    flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
    for s in iter(flattened_data):
        flattened_data[s] = [[data[s][x][y.slug] for y in numbering_schemes]+data[s][x]['seq'] for x in sorted(data[s])]
    
    context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
    context['segments'] = [x.slug for x in segments if len(data[x.slug])]
    context['data'] = flattened_data
    context['number_of_schemes'] = len(numbering_schemes)

    if download:
        raws = mutations.values('raw')
        rawmutations = MutationRaw.objects.filter(pk__in = raws).all()
        
        data = []
        for r in rawmutations:
            headers = []
            values = {}
            for field, val in r:
                headers.append(field)
                values[field] = val
            if values['id'] in mutations_display_generic_number:
                values['generic'] = mutations_display_generic_number[values['id']]
            else:
                values['generic'] = ''
            data.append(values)
        headers = ['reference', 'protein', 'mutation_pos', 'generic', 'mutation_from', 'mutation_to', 
        'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
        'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value', 
        'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_type', 'opt_wt',
        'opt_mu', 'opt_sign', 'opt_percentage', 'opt_qual','opt_agonist', 'added_date'                
         ] #'added_by',

        #CSV SOLUTION
        csv = "sep=;\n"
        for h in headers:
            csv += h + ";"
        csv += "\n"
        for d in data:
            for h in headers:
                csv += str(d[h]) + ";"
            csv += "\n"
        response_kwargs['content_type'] = 'text/csv'
        response = HttpResponse(**response_kwargs) 
        response['Content-Disposition'] = 'attachment; filename="GPCRdb_mutant_data.csv"'
        template = loader.get_template('mutation/csv.html')
        response.write(template.render(Context({'mutations': csv })))
        #return response

        #EXCEL SOLUTION
        output = BytesIO()
        workbook = xlsxwriter.Workbook(output)
        worksheet = workbook.add_worksheet()

        col = 0
        for h in headers:
            worksheet.write(0, col, h)
            col += 1
        row = 1
        for d in data:
            col = 0
            for h in headers:
                worksheet.write(row, col, d[h])
                col += 1
            row += 1
        workbook.close()
        output.seek(0)
        xlsx_data = output.read()

        response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename=GPCRdb_mutant_data.xlsx' #% 'mutations'
        return response


    else:        
        return render(request, 'mutation/list.html', {'mutations': mutations, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot, 'data':context['data'], 
            'header':context['header'], 'segments':context['segments'], 'mutations_pos_list' : json.dumps(mutations_pos_list), 'number_of_schemes':len(numbering_schemes), 'protein_ids':str(protein_ids)})

# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")


# Create your views here.
def ajax(request, slug, **response_kwargs):
    if '[' in slug:
        x = ast.literal_eval(urllib.parse.unquote(slug))
        #print(x)
        #h = HTMLParser.HTMLParser()
        #print(h.unescape(slug))
        mutations = MutationExperiment.objects.filter(protein__pk__in=x).order_by('residue__sequence_number').prefetch_related('residue')
    else:
        mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
    #print(mutations)
    #return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for mutation in mutations:
        if mutation.residue.sequence_number not in jsondata: jsondata[mutation.residue.sequence_number] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        jsondata[mutation.residue.sequence_number].append([mutation.foldchange,ligand,qual])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def ajaxSegments(request, slug, segments, **response_kwargs):
    print('hi')
    if '[' in slug:
        print("OMG OMG OMG")
        print(urllib.parse.unquote(slug))
        x = ast.literal_eval(urllib.parse.unquote(slug))
        segments = ast.literal_eval(urllib.parse.unquote(segments))
        print(segments)
        #h = HTMLParser.HTMLParser()
        #print(h.unescape(slug))
        mutations = MutationExperiment.objects.filter(protein__pk__in=x,residue__protein_segment__slug__in=segments).order_by('residue__sequence_number').prefetch_related('residue')
        print(mutations)
    else:
        mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
    #print(mutations)
    #return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for mutation in mutations:
        if mutation.residue.sequence_number not in jsondata: jsondata[mutation.residue.sequence_number] = []
        jsondata[mutation.residue.sequence_number].append(mutation.foldchange)

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def importmutation(request):

    rows = loaddatafromexcel('/vagrant/protwis/mutation/import.xlsx')

    rows = analyse_rows(rows)

    whattoreturn = []
    c = 0
    skipped = 0
    inserted = 0
    for r in rows:
        print(inserted,skipped)
        #print(r)
        raw_id = insert_raw(r)
        ref_id = check_reference(r['reference'])
        lig_id = get_ligand(r)

        if (ref_id==None):
            whattoreturn.append(['Skipped due to no working reference',r])
            skipped += 1
            continue

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

        check=Residue.objects.filter(protein_conformation__protein=protein_id,sequence_number=r['mutation_pos'])
        if check.exists():
            check=Residue.objects.get(protein_conformation__protein=protein_id,sequence_number=r['mutation_pos'])
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

        obj, created = MutationExperimentalType.objects.get_or_create(type=r['exp_type'])
        exp_type_id = obj

        obj, created = MutationFunc.objects.get_or_create(func=r['exp_func'])
        exp_func_id = obj

        obj, created = MutationMeasure.objects.get_or_create(measure=r['exp_mu_effect_type'])
        exp_measure_id = obj

        obj, created = MutationQual.objects.get_or_create(qual=r['exp_mu_effect_qual'], prop=r['exp_mu_effect_ligand_prop'])
        exp_qual_id = obj

        obj, created = MutationLigandRef.objects.get_or_create(reference=r['exp_mu_ligand_ref'])
        effect_ligand_reference_id = obj

        obj, created =  MutationOptional.objects.get_or_create(type=r['opt_type'], wt=r['opt_wt'], mu=r['opt_mu'], sign=r['opt_sign'], percentage=r['opt_percentage'], qual=r['opt_qual'], agonist=r['opt_agonist'])
        exp_opt_id = obj

        obj, created =  Mutation.objects.get_or_create(amino_acid=r['mutation_to'],protein=protein_id, residue=residue_id)
        mutation_id = obj


        
        logtypes = ['pEC50','pIC50','pK']
        
        
        foldchange = 0
        typefold = ''
        if r['exp_mu_effect_type']=='Activity/affinity' and r['exp_wt_value']!=0:
                    
            if re.match("(" + ")|(".join(logtypes) + ")", r['exp_type']):  #-log values!
                foldchange = round(math.pow(10,-r['exp_mu_value_raw'])/pow(10,-r['exp_wt_value']),3);
                typefold = r['exp_type']+"_log"
            else:
                foldchange = round(r['exp_mu_value_raw']/r['exp_wt_value'],3);
                typefold = r['exp_type']+"_not_log"
            
            
            if foldchange<1 and foldchange!=0:
                foldchange = -round((1/foldchange),3)
            elif r['exp_mu_effect_type'] =='Fold effect (mut/wt)':
                foldchange = round(r['exp_mu_value_raw'],3);
                if foldchange<1: foldchange = -round((1/foldchange),3);
        

        obj, created = MutationExperiment.objects.get_or_create(
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

        mutation=mutation_id, 
        wt_value=r['exp_wt_value'], #
        wt_unit=r['exp_wt_unit'], 

        mu_value = r['exp_mu_value_raw'],
        mu_sign = r['exp_mu_effect_sign'], 
        foldchange = foldchange
        #foldchange = 1

        #added_by='munk', 
        #added_date=datetime.now()
        )
        #print(foldchange)
        mut_id = obj.id


        #whattoreturn.append([protein_id,residue_id,raw_id,ref_id,lig_id,ligclass_id,exp_type_id,exp_func_id,exp_measure_id,exp_qual_id,typefold,foldchange,mut_id])
        whattoreturn.append(['Inserted',protein_id.entry_name,residue_id.sequence_number,foldchange])
        inserted += 1
        c += 1
        #if c>10: break

    print('DONE!')
    context = {'rows': whattoreturn, 'skipped' : skipped, 'inserted' : inserted}

    #return HttpResponse(ref_id)
    return render(request,'mutation/index.html',context)