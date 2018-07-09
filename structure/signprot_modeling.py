

class SignprotModeling():
	def __init__(self, main_structure, signprot, template_source):
		pass

	def run(self):
		
        self.signprot_complex = SignprotComplex.objects.get(structure=self.main_structure)
        structure_signprot= self.signprot_complex.protein
        if self.signprot!=False:
            self.target_signprot = Protein.objects.get(entry_name=self.signprot)
        else:
            self.target_signprot = self.signprot_complex.protein
        self.signprot_protconf = ProteinConformation.objects.get(protein=self.target_signprot)
        sign_a = GProteinAlignment()
        sign_a.run_alignment(self.target_signprot)
        io = StringIO(self.main_structure.pdb_data.pdb)
        assign_cgn = as_gn.GenericNumbering(pdb_file=io, pdb_code=self.main_structure.pdb_code.index, sequence_parser=True, signprot=structure_signprot)
        signprot_pdb_array = assign_cgn.assign_cgn_with_sequence_parser(self.signprot_complex.chain)
        new_array = OrderedDict()

        # Initiate complex part of template source
        source_resis = Residue.objects.filter(protein_conformation__protein=self.target_signprot)
        for res in source_resis:
            if res.protein_segment.slug not in self.template_source:
                self.template_source[res.protein_segment.slug] = OrderedDict()
            if res.protein_segment.category=='loop':
                self.template_source[res.protein_segment.slug][str(res.sequence_number)] = [None, None]
            else:
                self.template_source[res.protein_segment.slug][res.display_generic_number.label] = [self.main_structure, self.main_structure]

        # Superimpose missing regions H1 - hfs2
        alt_complex_struct = None
        segs_for_alt_complex_struct = []
        if self.main_structure.pdb_code.index!='3SN6':
            segs_for_alt_complex_struct = ['H1', 'h1ha', 'HA', 'hahb', 'HB', 'hbhc', 'HC', 'hchd', 'HD', 'hdhe', 'HE', 'hehf', 'HF', 'hfs2']
            alt_complex_struct = Structure.objects.get(pdb_code__index='3SN6')
            io = StringIO(alt_complex_struct.pdb_data.pdb)
            alt_signprot_complex = SignprotComplex.objects.get(structure__pdb_code__index='3SN6')
            alt_assign_cgn = as_gn.GenericNumbering(pdb_file=io, pdb_code='3SN6', sequence_parser=True, signprot=alt_signprot_complex.protein)
            alt_signprot_pdb_array = alt_assign_cgn.assign_cgn_with_sequence_parser(alt_signprot_complex.chain)
            before_cgns = ['G.HN.50', 'G.HN.51', 'G.HN.52', 'G.HN.53']
            after_cgns =  ['G.H5.03', 'G.H5.04', 'G.H5.05', 'G.H5.06']
            orig_residues1 = parse.fetch_residues_from_array(signprot_pdb_array['HN'], before_cgns)
            orig_residues2 = parse.fetch_residues_from_array(signprot_pdb_array['H5'], after_cgns)
            orig_residues = parse.add_two_ordereddict(orig_residues1, orig_residues2)

            alt_residues1 = parse.fetch_residues_from_array(alt_signprot_pdb_array['HN'], before_cgns)
            alt_residues2 = parse.fetch_residues_from_array(alt_signprot_pdb_array['H5'], after_cgns)
            alt_middle = OrderedDict()
            for s in segs_for_alt_complex_struct:
                alt_middle = parse.add_two_ordereddict(alt_middle, alt_signprot_pdb_array[s])
                self.update_template_source(list(self.template_source[s].keys()), alt_complex_struct, s)

            alt_residues = parse.add_two_ordereddict(parse.add_two_ordereddict(alt_residues1, alt_middle), alt_residues2)
            del_list = []
            for r, t in alt_middle.items():
                if t=='x':
                    del_list.append(r)
            for r in del_list:
                del alt_residues[r]
            
            superpose = sp.LoopSuperpose(orig_residues, alt_residues)
            new_residues = superpose.run()
            key_list = list(new_residues.keys())[4:-4]
            for key in key_list:
                seg = key.split('.')[1]
                signprot_pdb_array[seg][key] = new_residues[key]

            # Let Modeller model buffer regions
            trimmed_residues.append('s1h1_6')
            trimmed_residues.append('G.S2.01')
            trimmed_residues.append('G.S2.02')
            trimmed_residues.append('s4h3_4')
            trimmed_residues.append('s4h3_5')

        # New loop alignments for signprot. If length differs between ref and temp, buffer is created in the middle of the loop
        loops = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='Gprotein', category='loop')]
        loops_to_model = []
        for r_seg, t_seg, a_seg in zip(sign_a.reference_dict, sign_a.template_dict, sign_a.alignment_dict):
            if r_seg in loops:
                loop_length = len(sign_a.reference_dict[r_seg])
                ref_loop = [i for i in list(sign_a.reference_dict[r_seg].values()) if i not in ['x','-']]
                ref_keys = [i for i in list(sign_a.reference_dict[r_seg].keys()) if i not in ['x','-']]
                ref_loop_residues = Residue.objects.filter(protein_conformation__protein=self.target_signprot, protein_segment__slug=r_seg)
                temp_loop = [i for i in list(sign_a.template_dict[t_seg].values()) if i not in ['x','-']]
                temp_keys = [i for i in list(sign_a.template_dict[t_seg].keys()) if i not in ['x','-']]
                temp_loop_residues = Residue.objects.filter(protein_conformation__protein=structure_signprot, protein_segment__slug=r_seg)
                ref_out, temp_out, align_out = OrderedDict(), OrderedDict(), OrderedDict()
                print(r_seg)
                print(ref_loop)
                print(temp_loop)
                # ref is longer
                if len(ref_loop)>len(temp_loop):
                    mid_temp = math.ceil(len(temp_loop)/2)
                    j = 0
                    for i in range(0, loop_length):
                        key = r_seg+'_'+str(i+1)
                        if i+1<=mid_temp:
                            temp_out[key] = temp_loop[i]
                        elif mid_temp<i+1<=loop_length-mid_temp+1:
                            if i+1==loop_length-mid_temp+1 and len(temp_loop)%2==0:
                                temp_out[key] = temp_loop[mid_temp+j]


                                # self.compare_and_update_template_source(r_seg, segs_for_alt_complex_struct, signprot_pdb_array, temp_keys, mid_temp, j, ref_loop_residues, i, alt_complex_struct)


                                j+=1
                            else:
                                temp_out[key.replace('_','?')] = '-'
                                self.template_source[r_seg][str(ref_loop_residues[mid_temp+j].sequence_number)] = [None, None]
                        else:
                            temp_out[key] = temp_loop[mid_temp+j]
                            j+=1
                    for i, j in enumerate(list(sign_a.reference_dict[r_seg].values())):
                        key = r_seg+'_'+str(i+1)
                        try:
                            temp_out[key]
                            ref_out[key] = j
                        except:
                            ref_out[key.replace('_','?')] = j
                            # self.template_source[r_seg][str(ref_loop_residues[i].sequence_number)] = [None, None]
                        i+=1
                # temp is longer
                elif len(ref_loop)<len(temp_loop):
                    mid_ref = math.ceil(len(ref_loop)/2)
                    j = 0
                    for i in range(0, loop_length):
                        key = r_seg+'_'+str(i+1)
                        print(key, i, j, mid_ref)
                        if i+1<=mid_ref:
                            ref_out[key] = ref_loop[i]
                            # if r_seg in segs_for_alt_complex_struct and signprot_pdb_array[r_seg][:
                            #     self.template_source[r_seg][str(ref_loop_residues[i].sequence_number)] = [alt_complex_struct, alt_complex_struct]
                            # else:
                            #     self.template_source[r_seg][str(ref_loop_residues[i].sequence_number)] = [self.main_structure, self.main_structure]
                        elif mid_ref<i+1<=loop_length-mid_ref+1:
                            if i+1==loop_length-mid_ref+1 and len(ref_loop)%2==0:
                                ref_out[key] = ref_loop[mid_ref+j]

                                # self.compare_and_update_template_source(r_seg, segs_for_alt_complex_struct, signprot_pdb_array, ref_keys, mid_ref, j, temp_loop_residues, i, alt_complex_struct)

                                j+=1
                            else:
                                ref_out[key.replace('_','?')] = '-'
                                self.template_source[r_seg][str(ref_loop_residues[mid_ref+j].sequence_number)] = [None, None]
                        else:
                            ref_out[key] = ref_loop[mid_ref+j]
                            # if r_seg in segs_for_alt_complex_struct:
                            #     self.template_source[r_seg][str(ref_loop_residues[mid_ref+j].sequence_number)] = [alt_complex_struct, alt_complex_struct]
                            # else:
                            #     self.template_source[r_seg][str(ref_loop_residues[mid_ref+j].sequence_number)] = [self.main_structure, self.main_structure]
                            j+=1
                    for i, j in enumerate(list(sign_a.template_dict[t_seg].values())):
                        key = r_seg+'_'+str(i+1)
                        try:
                            ref_out[key]
                            temp_out[key] = j
                        except:
                            temp_out[key.replace('_','?')] = j
                            # self.template_source[r_seg][str(ref_loop_residues[i].sequence_number)] = [None, None]
                        i+=1
                    loops_to_model.append(r_seg)
                # ref and temp length equal
                else:
                    c = 1
                    for i, j in zip(list(sign_a.reference_dict[r_seg].values()), list(sign_a.template_dict[t_seg].values())):
                        ref_out[r_seg+'_'+str(c)] = i
                        temp_out[r_seg+'_'+str(c)] = j
                        # if r_seg in segs_for_alt_complex_struct:
                        #     self.template_source[r_seg][str(ref_loop_residues[c-1].sequence_number)] = [alt_complex_struct, alt_complex_struct]                            
                        # else:
                        #     self.template_source[r_seg][str(ref_loop_residues[c-1].sequence_number)] = [self.main_structure, self.main_structure]
                        c+=1
                c = 1
                # update alignment dict
                for i, j in zip(list(ref_out.values()), list(temp_out.values())):
                    key = r_seg+'_'+str(c)
                    if i=='-' or j=='-':
                        align_out[key.replace('_','?')] = '-'
                    elif i!=j:
                        align_out[key] = '.'
                    elif i==j:
                        align_out[key] = i
                    c+=1

                # update pdb array
                new_pdb_array = OrderedDict()
                atoms_list = list(signprot_pdb_array[t_seg].values())
                j = 0
                for t_c, t in temp_out.items():
                    if t!='-':
                        for i in range(j, len(atoms_list)):
                            if atoms_list[j]!='-':
                                new_pdb_array[t_c] = atoms_list[j]
                                break
                        j+=1
                    else:
                        new_pdb_array[t_c] = 'x'
                print(r_seg)
                print(ref_out)
                print(temp_out)
                print(align_out)
                # update dictionary keys with '?' if no backbone template
                ref_out_final, temp_out_final, align_out_final, new_pdb_array_final = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
                # self.template_source[r_seg] = OrderedDict()
                for i,j in new_pdb_array.items():
                    if '?' not in i and j=='x':
                        ref_out_final[i.replace('_','?').replace('.','?')] = ref_out[i]
                        temp_out_final[i.replace('_','?').replace('.','?')] = temp_out[i]
                        align_out_final[i.replace('_','?').replace('.','?')] = align_out[i]
                        new_pdb_array_final[i.replace('_','?').replace('.','?')] = new_pdb_array[i]
                        # self.template_source[r_seg][i.replace('_','?').replace('.','?')] = [None, None]
                    else:
                        ref_out_final[i] = ref_out[i]
                        temp_out_final[i] = temp_out[i]
                        align_out_final[i] = align_out[i]
                        new_pdb_array_final[i] = new_pdb_array[i]
                        # if r_seg not in segs_for_alt_complex_struct:
                        #     self.template_source[r_seg][i] = [self.main_structure, self.main_structure]
                        # else:
                        #     self.template_source[r_seg][i] = [alt_complex_struct, alt_complex_struct]


                sign_a.reference_dict[r_seg] = ref_out_final
                sign_a.template_dict[t_seg] = temp_out_final
                sign_a.alignment_dict[a_seg] = align_out_final
                signprot_pdb_array[r_seg] = new_pdb_array_final
                
                align_loop = list(sign_a.alignment_dict[a_seg].values())
        
        # for i,j,k,l in zip(sign_a.reference_dict, sign_a.template_dict, sign_a.alignment_dict, signprot_pdb_array):
        #     for v,b,n,m in zip(sign_a.reference_dict[i], sign_a.template_dict[j], sign_a.alignment_dict[k], signprot_pdb_array[l]):
        #         print(v, b, n, m, sign_a.reference_dict[i][v], sign_a.template_dict[j][b], sign_a.alignment_dict[k][n], signprot_pdb_array[l][m])


        for seg, values in sign_a.reference_dict.items():
            new_array[seg] = OrderedDict()
            # self.template_source[seg] = OrderedDict()
            for key, res in values.items():
                try:
                    if signprot_pdb_array[seg][key]=='x':
                        new_array[seg][key] = 'x'
                        try:
                            self.template_source[seg][key] = [None, None]
                        except:
                            pass
                    else:
                        new_array[seg][key] = signprot_pdb_array[seg][key]
                except:
                    if res!='-':
                        new_array[seg][key] = '-'
                        try:
                            self.template_source[seg][key] = [None, None]
                        except:
                            pass
            a.reference_dict[seg] = values
        for seg, values in sign_a.template_dict.items():
            for key, res in values.items():
                if new_array[seg][key]=='x':
                    values[key] = 'x'
                else:
                    pdb_res = PDB.Polypeptide.three_to_one(new_array[seg][key][0].get_parent().get_resname())
                    if pdb_res!=values[key]:
                        values[key] = pdb_res
            a.template_dict[seg] = values
        for seg, values in sign_a.alignment_dict.items():
            for key, res in values.items():
                if new_array[seg][key]=='x':
                    values[key] = 'x'
            a.alignment_dict[seg] = values
        signprot_pdb_array = new_array
        for seg, values in signprot_pdb_array.items():
            main_pdb_array[seg] = values

        delete_HN_begin = []
        for i in a.reference_dict['HN']:
            if i=='G.HN.30':
                break
            delete_HN_begin.append(i)
        for d in delete_HN_begin:
            del a.reference_dict['HN'][d]
            del a.template_dict['HN'][d]
            del a.alignment_dict['HN'][d]
            del main_pdb_array['HN'][d]
            try:
                del self.template_source['HN'][d]
            except:
                pass

        # add residues to model to trimmed_residues
        gprot_segments = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='Gprotein')]
        for i,j in a.reference_dict.items():
            if i in gprot_segments:
                for k,l in j.items():
                    if '?' in k or main_pdb_array[i][k] in ['-','x']:
                        trimmed_residues.append(k)
                    if i in loops_to_model:
                        trimmed_residues.append(k)

        # custom mods
        long_HG_prots = Protein.objects.filter(family__name='Gs')
        if structure_signprot in long_HG_prots and self.target_signprot not in long_HG_prots:
            trimmed_residues.append('G.HG.08')
            trimmed_residues.append('G.HG.09')
            trimmed_residues.append('G.HG.12')
            trimmed_residues.append('G.HG.13')
            trimmed_residues.append('G.HG.14')
            trimmed_residues.append('G.HG.16')
            trimmed_residues.append('G.HG.17')
        if 'hgh4' in loops_to_model:
            trimmed_residues.append('G.H4.01')
            trimmed_residues.append('G.H4.02')
            trimmed_residues.append('G.H4.03')

	    # for i,j,k,l in zip(sign_a.reference_dict, sign_a.template_dict, sign_a.alignment_dict, signprot_pdb_array):
	    #     pprint.pprint(self.template_source[i])
	    #     for v,b,n,m in zip(sign_a.reference_dict[i], sign_a.template_dict[j], sign_a.alignment_dict[k], signprot_pdb_array[l]):
	    #         print(v, b, n, m, sign_a.reference_dict[i][v], sign_a.template_dict[j][b], sign_a.alignment_dict[k][n], signprot_pdb_array[l][m])