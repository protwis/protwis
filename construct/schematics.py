from protein.models import Protein, ProteinConformation
from residue.models import Residue
from structure.models import Structure

import time

def generate_schematic(c):
    global fusion_protein_name
    # print("GENERATING!")

    time1 = time.time()

    ### PREPARE DATA
    summary = {}
    annotations = {}
    json_annotations = {}

    fusion_position,fusion_result, linkers = c.fusion()


    time2 = time.time()
    # print('%s function took %0.3f ms' % ('fusion', (time2-time1)*1000.0))

    if fusion_result:
        fusion_protein_name = fusion_result[0][2]
    else:
        fusion_protein_name = None

    summary['solubilization'] = ''
    if c.solubilization:
        for chem in c.solubilization.chemical_list.chemicals.prefetch_related('chemical').all():
            summary['solubilization'] += """{} ({} {})<br>""".format(chem.chemical.name,chem.concentration,chem.concentration_unit)
    time2 = time.time()
    # print('%s function took %0.3f ms' % ('summaries', (time2-time1)*1000.0))


    summary['purification'] = ''
    if c.purification:
        for step in c.purification.steps.all():
            if step.name!='None':
                summary['purification'] += """{}<br>""".format(step.name)
    time2 = time.time()
    # print('%s function took %0.3f ms' % ('summaries', (time2-time1)*1000.0))


    summary['crystallization_chems'] = ''
    if c.crystallization:
        for clist in c.crystallization.chemical_lists.all():
            # summary['crystallization_chems'] += """{}<br>""".format(clist.name)    
            for chem in clist.chemicals.prefetch_related('chemical').all():
                summary['crystallization_chems'] += """{} ({} {})<br>""".format(chem.chemical.name,chem.concentration,chem.concentration_unit)


    time2 = time.time()
    # print('%s function took %0.3f ms' % ('summaries', (time2-time1)*1000.0))

    n_term = {}
    c_term = {}
    insert = {}
    deletion = {}
    pair_no = 0
    summary['modifications'] = ''
    for mod in c.modifications.all():
        if mod.position_type == 'pair':
            pair_no += 1
            pair_info = "Pair No"+str(pair_no)+" - "+str(mod.pos_start)+" : "+str(mod.pos_end)
            pair_short = str(pair_no)
        else:
            pair_info = str(mod.pos_start)
            pair_short = ''

        if mod.modification=='Disulfide bond': continue

        annotations[mod.pos_start] = ['mod','modification','Modified<br>'+mod.modification+"<br>"+pair_info,mod,pair_short]
        annotations[mod.pos_end] = ['mod','modification','Modified<br>'+mod.modification+"<br>"+pair_info,mod,pair_short]
        json_annotations[mod.pos_start] = ['mod','Modified<br>'+mod.modification+"<br>"+pair_info,"yellow","black"]
        json_annotations[mod.pos_end] = ['mod','Modified<br>'+mod.modification+"<br>"+pair_info,"yellow","black"]
        summary['modifications'] += pair_info+" - "+mod.modification+"<br>"


    summary['mutations'] = ''
    for mut in c.mutations.all():
        annotations[mut.sequence_number] = ['mut','mutation','Mutated from '+mut.wild_type_amino_acid+' to '+mut.mutated_amino_acid,mut]
        json_annotations[mut.sequence_number] = ['mut','Mutated from '+mut.wild_type_amino_acid+' to '+mut.mutated_amino_acid,"gold","black"]
        summary['mutations'] += mut.wild_type_amino_acid+str(mut.sequence_number)+mut.mutated_amino_acid+'<br>'

    summary['insertions'] = '<div style="text-align:left">'
    for aux in c.insertions.order_by('position').prefetch_related('insert_type').all():
        position_without_number = aux.position.split("_")[0]
        if aux.start:
            if aux.start==19 and aux.insert_type.name=='auto':
                continue
            if aux.start in insert:    
                print("ERROR Multiple inserts at same position",aux.start,c.name,aux)
                # if insert[aux.start].insert_type.name =='fusion':
                #     #continue to next
                #     # continue
                #     pass
            else:
                insert[aux.start] = []
            for i in range(aux.start,aux.end+1):
                annotations[i] = [aux.insert_type.name,'Insertion<br>Protein_type: '+aux.insert_type.name,aux]
                json_annotations[i] = ['ins','Insertion<br>Protein_type: '+aux.insert_type.name,"purple","white"]
            
            insert[aux.start].append(aux)

        if aux.position.startswith('N-term'):
            n_term[aux.position] = aux
        if aux.position.startswith('C-term'):
            c_term[aux.position] = aux
        summary['insertions'] += """{} - {} ({})<br>""".format(position_without_number,aux.insert_type.subtype,aux.insert_type.name)
    summary['insertions'] += '</div>'
    summary['deletions'] = ''
    for dele in c.deletions.all():
        deletion[dele.start] = ['del','deletion','Deleted from '+str(dele.start)+' to '+str(dele.end),dele]
        deletion[dele.end] = ['del','deletion','Deleted from '+str(dele.start)+' to '+str(dele.end),dele]

        summary['deletions'] += str(dele.start)+ " to "+str(dele.end)+'<br>'
        for i in range(dele.start,dele.end+1):
            annotations[i] = ['del','deletion','Deleted from '+str(dele.start)+' to '+str(dele.end),dele]
            json_annotations[i] = ['del','Deleted from '+str(dele.start)+' to '+str(dele.end),"red","white"]


    time2 = time.time()
    # print('%s function took %0.3f ms' % ('2', (time2-time1)*1000.0))

    # get residues
    residues = Residue.objects.filter(protein_conformation__protein=c.protein).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')


    time2 = time.time()
    # print('%s function took %0.3f ms' % ('3', (time2-time1)*1000.0))

    residues_custom = []
    residues_lookup = {}
    for r in residues:
        residues_lookup[r.sequence_number] = r
        r.print_amino_acid = r.amino_acid
        if r.sequence_number in annotations:
            if annotations[r.sequence_number][1] == 'deletion':
                continue
            if annotations[r.sequence_number][1] == 'insertion':
                continue
            if annotations[r.sequence_number][1] == 'mutation':
                r.print_amino_acid = annotations[r.sequence_number][3].mutated_amino_acid
        residues_custom.append(r)

    time2 = time.time()
    # print('%s function took %0.3f ms' % ('4', (time2-time1)*1000.0))

    results = {}
    results['annotations'] = json_annotations
    #### VERSION 1 SCHEMATIC

    chunk_size = 10
    r_chunks_schematic = []
    r_chunks_schematic_construct = []
    r_buffer = []
    a_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    ii = 0
    # create schemtics with annotations
    for i, r in enumerate(residues):
        if (r.protein_segment.slug != last_segment and i!=0) or i == len(residues)-1:
            a_list = {}
            if ii<20:
                counter = ii
                width = round(100/ii,1)
            else:
                counter = 20
                width = 5
            for a in a_buffer:
                temp = a[0]*counter // (ii)
                a_list[temp]  = a[1]
            if (last_segment):
                for a in range(counter):
                    if a<len(prev_r.protein_segment.slug):
                        title_cell_skip = 1
                    else:
                        title_cell_skip = 0

                    if a in a_list:
                        annotation = a_list[a]
                    else:
                        annotation = ''

                    if a==0:
                        r_buffer.append([prev_r, True, 0,annotation,width])
                    else:
                        r_buffer.append([prev_r, False, 0,annotation,width])

                r_chunks_schematic.append(r_buffer)
                #print(last_segment,prev_r.protein_segment.slug,counter,len(r_buffer),width,len(r_buffer)*width)
            border = True
            r_buffer = []
            a_buffer = []
            ii = 0

        if r.sequence_number in annotations:
            a_buffer.append([ii,annotations[r.sequence_number]])
        prev_r = r
        last_segment = r.protein_segment.slug
        ii+=1

    results['schematic_1_wt'] = r_chunks_schematic

    #### VERSION 1 SCHEMATIC CONSTRUC

    time2 = time.time()
    # print('%s function took %0.3f ms' % ('1', (time2-time1)*1000.0))

    for aux in sorted(n_term):
        name = n_term[aux].insert_type.subtype[:6]
        title_cell_skip = 0
        for i in range(20):
            if i==0:
                r_buffer.append([None, name, title_cell_skip,['insert','insertion',n_term[aux]]])
            else:
                if i<len(name):
                    title_cell_skip = 1
                else:
                    title_cell_skip = 0
                r_buffer.append([None, False, title_cell_skip,['insert','insertion',n_term[aux]]])
        r_chunks_schematic_construct.append(r_buffer)
        r_buffer = []

        ii = 0

    last_segment = False
    a_buffer = []
    r_buffer = []
    nudge = 0
    prev_r = None
    # create schemtics with annotations
    for i, r in enumerate(residues_custom):
        # print(i)
        if r.sequence_number-1 in deletion and i==0:
            r_buffer.append([None, False, 0,annotations[r.sequence_number-1]])
            nudge = 1

        if r.sequence_number+1 in deletion:
            if i < len(residues_custom)-1: #not last
                a_buffer.append([ii+1,annotations[r.sequence_number+1]])
            else:
                a_buffer.append([ii,annotations[r.sequence_number+1]])

        if (r.protein_segment.slug != last_segment and i!=0) or i == len(residues_custom)-1:
            a_list = {}
            for a in a_buffer:
                temp = a[0]*20 // (ii)
                a_list[temp]  = a[1]
            if (last_segment):
                for a in range(21):
                    if (a-nudge)<len(prev_r.protein_segment.slug):
                        title_cell_skip = 1
                    else:
                        title_cell_skip = 0

                    if a in a_list:
                        annotation = a_list[a]
                        if annotation[1] == 'deletion' and a==20:
                            no_r = True
                        else:
                            no_r = False
                    else:
                        annotation = ''
                        no_r = False

                    if a==0:
                        r_buffer.append([prev_r, prev_r.protein_segment.slug, 0,annotation])
                    else:
                        if no_r:
                            r_buffer.append([None, False, title_cell_skip,annotation])
                        else:
                            r_buffer.append([prev_r, False, title_cell_skip,annotation])
                r_chunks_schematic_construct.append(r_buffer)
                # print("2",last_segment,prev_r.protein_segment.slug,counter,len(r_buffer),width,len(r_buffer)*width)
            last_segment = r.protein_segment.slug
            border = True
            r_buffer = []
            a_buffer = []
            ii = 0
            nudge = 0
        #print(r_chunks_schematic_construct)

        if r.sequence_number in annotations:
            a_buffer.append([ii,annotations[r.sequence_number]])

        if prev_r:
            if i+1==len(residues_custom):
                prev_r = r
            if prev_r.sequence_number+1 in insert:
                for temp in insert[prev_r.sequence_number+1]:
                    name = temp.insert_type.subtype[:6]
                    title_cell_skip = 0
                    for i in range(21):
                        if i==0:
                            r_buffer.append([None, name, title_cell_skip,['insert','insertion',temp]])
                        else:
                            if i<len(name):
                                title_cell_skip = 1
                            else:
                                title_cell_skip = 0
                            r_buffer.append([None, False, title_cell_skip,['insert','insertion',temp]])
                    r_chunks_schematic_construct.append(r_buffer)
                    r_buffer = []


        if i+1==len(residues_custom) and r.sequence_number+1 in insert:
            name = insert[r.sequence_number+1].insert_type.subtype[:6]
            title_cell_skip = 0
            for i in range(21):
                if i==0:
                    #print(last_segment,name,len(r_chunks_schematic_construct))
                    r_buffer.append([None, name, title_cell_skip,['insert','insertion',insert[r.sequence_number+1]]])
                else:
                    if i<len(name):
                        title_cell_skip = 1
                    else:
                        title_cell_skip = 0
                    r_buffer.append([None, False, title_cell_skip,['insert','insertion',insert[r.sequence_number+1]]])
            r_chunks_schematic_construct.append(r_buffer)
            r_buffer = []


        last_segment = r.protein_segment.slug
        prev_r = r
        ii+=1

    for aux in sorted(c_term):
        # print("new",aux,c_term[aux])
        name = c_term[aux].insert_type.subtype[:6]
        title_cell_skip = 0
        for i in range(20):
            if i==0:
                r_buffer.append([None, name, title_cell_skip,['insert','insertion',c_term[aux]]])
            else:
                if i<len(name):
                    title_cell_skip = 1
                else:
                    title_cell_skip = 0
                r_buffer.append([None, False, title_cell_skip,['insert','insertion',c_term[aux]]])
        r_chunks_schematic_construct.append(r_buffer)
        r_buffer = []

        ii = 0


    results['schematic_1_c'] = r_chunks_schematic_construct
    # print(r_chunks_schematic_construct)

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False
        
        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True
        
        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)
        
        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []
        
        if r.sequence_number in annotations:
            annotation = annotations[r.sequence_number]
        else:
            annotation = None

        r_buffer.append((r, segment_title, title_cell_skip,annotation))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)


    results['residues_wt'] = r_chunks

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    r_chunks_custom = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0

    nudge = 0

    fusion = False

    for aux in sorted(n_term):
        name = n_term[aux].insert_type.subtype[:11]
        title_cell_skip = 0
        for i in range(chunk_size):
            if i==0:
                r_buffer.append([None, name, title_cell_skip,['insert',n_term[aux].insert_type.name,n_term[aux]]])
            else:
                if i<len(name):
                    title_cell_skip = 1
                else:
                    title_cell_skip = 0
                r_buffer.append([None, False, title_cell_skip,['insert',n_term[aux].insert_type.name,n_term[aux]]])
        r_chunks_custom.append(r_buffer)
        r_buffer = []

    title_cell_skip = 0
    for i, r in enumerate(residues_custom):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        if r.sequence_number-1 in deletion and i==0:
            # r_buffer.append((r, segment_title, title_cell_skip,deletion[r.sequence_number+1]))
            r_buffer.append((None, segment_title, title_cell_skip,deletion[r.sequence_number-1]))
            nudge -= 1
     
        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True
        
        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - (i-nudge) % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)
        
        if len(r_buffer) and i and (i-nudge) % chunk_size == 0:
            r_chunks_custom.append(r_buffer)
            r_buffer = []
        
        if r.sequence_number in annotations:
            annotation = annotations[r.sequence_number]
        else:
            annotation = None
        r_buffer.append((r, segment_title, title_cell_skip,annotation))

        if r.sequence_number+1 in deletion:
            # r_buffer.append((r, segment_title, title_cell_skip,deletion[r.sequence_number+1]))
            r_buffer.append((None, segment_title, title_cell_skip,deletion[r.sequence_number+1]))
            if title_cell_skip > 0: title_cell_skip -=1
            nudge -= 1

        if r.sequence_number+1 in insert:
            for temp in insert[r.sequence_number+1]:
                if len(r_buffer):
                    r_chunks_custom.append(r_buffer)
                r_buffer = []
                # for _ in range(chunk_size):
                #     r_buffer.append([r, segment_title, title_cell_skip,['insert','insertion']])
                for a in range(chunk_size):
                    name = temp.insert_type.subtype[:11]
                    if a==0:
                        r_buffer.append([r, name, title_cell_skip,['insert','insertion',temp]])
                    else:
                        if a<len(name):
                            title_cell_skip = 1
                        else:
                            title_cell_skip = 0
                        r_buffer.append([r, False, title_cell_skip,['insert','insertion',temp]])
                r_chunks_custom.append(r_buffer)
                r_buffer = []
                nudge = ((i+nudge) % chunk_size)+1

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks_custom.append(r_buffer)

    r_buffer = []
    for aux in sorted(c_term):
        # print('adding cterm',aux,c_term[aux].insert_type)
        name = c_term[aux].insert_type.subtype[:11]
        title_cell_skip = 0
        for i in range(chunk_size):
            if i==0:
                r_buffer.append([None, name, title_cell_skip,['insert',c_term[aux].insert_type.name,c_term[aux]]])
            else:
                if i<len(name):
                    title_cell_skip = 1
                else:
                    title_cell_skip = 0
                r_buffer.append([None, False, title_cell_skip,['insert',c_term[aux].insert_type.name,c_term[aux]]])
        r_chunks_custom.append(r_buffer)
        r_buffer = []

    results['residues_c'] = r_chunks_custom

    #BUILDING


    ## SCHEMATIC WT

    wt_schematic = """<div class="row no-wrap" id="schematic_seq_wt">"""
    for rs in r_chunks_schematic:
        wt_schematic += """<div class="construct_div" style="">
                <table style='table-layout: fixed;width:50px' ><tr style="height:10px;">"""
        for r in rs:
            if r[3]:
                wt_schematic += """<td class="seqv seqv-segment no-wrap {}">
                            <div data-toggle="tooltip" data-placement="top" data-html="true"
                            title="{}">&nbsp;</div></td>""".format(r[3][1],r[3][2])
            else:
                wt_schematic += """<td class="seqv seqv-segment no-wrap">
                            &nbsp;
                            </td>"""
        wt_schematic += """</tr><tr>"""
        for r in rs:
            if not r[1] and r[2]:
                pass
            else:
                if r[0].protein_segment.id % 2 == 0:
                    extra = "bg-success"
                else:
                    extra = "bg-info"

                wt_schematic += """<td class="seqv seqv-segment no-wrap {}">""".format(extra)
                if r[1]:
                    wt_schematic += r[0].protein_segment.slug
                else:
                    wt_schematic += "&nbsp;"
                wt_schematic += "</td>"

        wt_schematic += "</tr></table></div>"
    wt_schematic += "</div>"


    ## SCHEMATIC CONSTRUCT

    c_schematic = """<div class="row no-wrap" id="schematic_seq_c">"""
    for rs in r_chunks_schematic_construct:
        c_schematic += """<div class="construct_div" style="">
                <table style='table-layout: fixed;width:50px'><tr style="height:10px;">"""
        for r in rs:
            if r[3]:
                c_schematic += """<td class="seqv seqv-segment no-wrap {}">
                            <div data-toggle="tooltip" data-placement="top" data-html="true"
                            title="{}">&nbsp;</div></td>""".format(r[3][1],r[3][2])
            else:
                c_schematic += """<td class="seqv seqv-segment no-wrap">
                            &nbsp;
                            </td>"""
        c_schematic += """</tr><tr>"""
        for r in rs:
            if not r[1] and r[2]:
                pass
            else:
                try:
                    extra = "bg-"+r[3][1]
                except:
                    extra = ""
                if r[0]: 
                    if r[0].protein_segment.id % 2 == 0:
                        extra = "bg-success"
                    else:
                        extra = "bg-info"

                colspan = 1

                if r[3] and r[1]:
                    text = r[1]
                    colspan = len(text)
                elif r[1]:
                    text = r[0].protein_segment.slug
                    colspan = len(text)
                else:
                    text = "&nbsp;"

                c_schematic += """<td colspan={} class="seqv seqv-segment no-wrap {}">{}</td>""".format(colspan,extra,text)

        c_schematic += "</tr></table></div>"
    c_schematic += "</div>"




    ### BUILD WT SPECIFIC TABLE
    wt_schematic_table = ""
    order_list = ['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','C-term'] #,'ICL4'
    i = 0
    for block in order_list:
        wt_schematic_table += "<td>"
        if i<len(r_chunks_schematic):
            slug = r_chunks_schematic[i][0][0].protein_segment.slug
            if block==slug:
                wt_schematic_table += create_block(r_chunks_schematic[i])
            else:
                i -=1 #not found match, dont move up
            i += 1
        wt_schematic_table += "</td>"


    results['schematic_2_wt'] = wt_schematic_table


    ### BUILD CONSTRUCT SPECIFIC TABLE
    c_schematic_table = ""
    order_list = ['pre','N-term','TM1','ICL1','TM2','ECL1','TM3','insert','ICL2','insert','TM4','ECL2','TM5','insert','ICL3','insert','TM6','ECL3','TM7','H8','C-term','post'] #,'ICL4'
    i = 0
    problem = False
    tm_corrected = False
    # print(r_chunks_schematic_construct)
    for block in order_list:
        c_schematic_table += "<td>"
        if problem:
            i -= problem
            problem = False
            tm_corrected = True
        # print(block,r_chunks_schematic_construct[i][0][1])
        try:
            if block=='TM5' and r_chunks_schematic_construct[i][0][1]!='TM5':
                # print('problem!',block,r_chunks_schematic_construct[i][0][1])
                problem = 1
                a = True
                while a and i<len(r_chunks_schematic_construct):
                    if r_chunks_schematic_construct[i][0][1]!='TM5':
                        i+=1
                        problem += 1
                    else:
                        a = False
                # print(r_chunks_schematic_construct[i][0][1],problem)
            if block=='ICL3' and tm_corrected and r_chunks_schematic_construct[i][0][1]!='ICL3':
                # print('try to find icl3')
                a = True
                while a and i<len(r_chunks_schematic_construct):
                    i += 1
                    if r_chunks_schematic_construct[i][0][1]=='ICL3':
                        # print('found ICL3',i)
                        a = False
            if block=='TM6' and tm_corrected and r_chunks_schematic_construct[i][0][1]!='TM6':
                # print('try to find TM6')
                a = True
                while a and i<len(r_chunks_schematic_construct):
                    i += 1
                    if r_chunks_schematic_construct[i][0][1]=='TM6':
                        # print('found TM6',i)
                        a = False
        except:
            pass
        if i<len(r_chunks_schematic_construct):
            if block=='pre':
                a = True
                c_schematic_table += "<div style2='float:right'><table align='right' width2='100%' class='no-wrap'><tr>"
                while a and i<len(r_chunks_schematic_construct):
                    if r_chunks_schematic_construct[i][0][3]:
                        if r_chunks_schematic_construct[i][0][3][0]=='insert':
                            c_schematic_table += "<td style='min-width:80px'>"+create_block(r_chunks_schematic_construct[i]) +"</td>"
                            i+=1
                        else:
                            a = False
                    else:
                        a = False
                c_schematic_table += "</tr></table></div>"
            elif block=='post':
                a = True
                c_schematic_table += "<div style2='float:left'><table align='left' width2='100%' class='no-wrap'><tr>"
                while a and i<len(r_chunks_schematic_construct):
                    if r_chunks_schematic_construct[i][0][3]:
                        if r_chunks_schematic_construct[i][0][3][0]=='insert' :
                            c_schematic_table += "<td style='min-width:80px'>"+create_block(r_chunks_schematic_construct[i]) +"</td>"
                            i+=1
                        else:
                            a = False
                    else:
                        a = False
                c_schematic_table += "</tr></table></div>"
            elif block=='insert':
                # print("insert?",r_chunks_schematic_construct[i])
                a = True
                c_schematic_table += "<div style2='float:right'><table align='right' width2='100%' class='no-wrap'><tr>"
                while a and i<len(r_chunks_schematic_construct):
                    if r_chunks_schematic_construct[i][0][3]:
                        if r_chunks_schematic_construct[i][0][3][0]=='insert' :
                            c_schematic_table += "<td style='min-width:80px'>"+create_block(r_chunks_schematic_construct[i]) +"</td>"
                            i+=1
                        else:
                            a = False
                    else:
                        a = False
                c_schematic_table += "</tr></table></div>"
            else:
                if block==r_chunks_schematic_construct[i][0][1]:
                    c_schematic_table += create_block(r_chunks_schematic_construct[i])
                    # print('found!',r_chunks_schematic_construct[i][0][1])
                elif block==r_chunks_schematic_construct[i][1][1]: #if started with deletion
                    c_schematic_table += create_block(r_chunks_schematic_construct[i])
                else:
                    # print('not found!',i,block,r_chunks_schematic_construct[i][0][1],r_chunks_schematic_construct[i][1][1])
                    i -=1 #not found match, dont move up
                    c_schematic_table += "&nbsp;"
                i += 1
        c_schematic_table += "</td>"


    results['schematic_2_c'] = c_schematic_table
        
    results['summary'] = summary

    time2 = time.time()
    # print('%s function took %0.3f ms' % ('final', (time2-time1)*1000.0))
    #print('done')
    results['schematic_1_c'] = '' #FAILS NOT REQUIRED THO
    results['schematic_1_wt'] = '' #FAILS NOT REQUIRED THO
    # results['residues_c'] = '' #FAILS
    # results['residues_wt'] = '' ## FAILS
    # results['schematic_2_wt'] = '' NOT FAILS
    # results['schematic_2_c'] = '' # Not fails
    # results['summary'] = '' # Not fails
    # results['annotations'] = ''# Not fails
    return results


def create_block(chunk):

    temp = """<div class="construct_div" style="">
            <table class='schematic-block'><tr>"""
    blank = 0
    blank_max = 20
    min_colspan = 1
    if len(chunk)<10:
        min_colspan = round(20/len(chunk))
    for i,r in enumerate(chunk):
        if r[3]:
            extra = r[3][1]
            if r[3][1]=='insertion' and i>19:
                continue #pass to prevent too many

            if r[3][1]=='insertion':
                r[3][1] = r[3][2].autotype()

                if r[3][2].insert_type.subtype==fusion_protein_name:
                    r[3][1] = 'fusion'

                extra = r[3][1]
                extra += "top"

            if blank: 
                temp += """<td class="seqv seqv-segment no-wrap" colspan="{}">&nbsp;</td>""".format(str(blank))
                blank = 0
            temp += """<td class="seqv seqv-segment no-wrap {}" colspan="{}">
                        <div data-toggle="tooltip" data-placement="top" data-html="true"
                        title="{}">&nbsp;</div></td>""".format(extra,min_colspan,r[3][2])
            blank_max = blank_max-i
        else:
            blank += 1
            if blank > blank_max: 
                blank = blank_max

    if blank: 
        temp += """<td class="seqv seqv-segment no-wrap" colspan="{}">&nbsp;</td>""".format(str(blank))
    temp += """</tr><tr>"""
    for r in chunk:
        if not r[1] and r[2] and 2==3:
            pass
        else:
            try:
                extra = "bg-"+r[3][1]
            except:
                extra = ""
            if r[0]: 
                if r[0].protein_segment.id % 2 == 0:
                    extra = "bg-success"
                else:
                    extra = "bg-info"

            if r[3] and r[1]:
                text = r[1]
                if text == True:
                    text = r[0].protein_segment.slug
            elif r[1]:
                text = r[0].protein_segment.slug
            else:
                text = "&nbsp;"
            if r[1] and r[3]:
                if r[3][0] == 'insert':
                    # print(text,r)
                    name = r[3][2].insert_type.subtype 
                    # print(name)
                    extra = r[3][1]
                    if name==fusion_protein_name:
                        extra = 'fusion'

            if r[1]:
                temp += """<td class="seqv seqv-segment no-wrap {} " colspan=20>{}</td>""".format(extra,text)
            elif r[0]==None and r[3] and r[3][0]!='insert':
                temp += """<td class="seqv seqv-segment no-wrap {}" colspan=1>{}</td>""".format(extra,text)

    temp += "</tr></table></div>"
    return temp