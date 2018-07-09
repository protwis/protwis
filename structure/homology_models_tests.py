from structure.models import Structure

from collections import OrderedDict
import Bio.PDB as PDB


class HomologyModelsTests(object):
    
    def pdb_alignment_mismatch(self, alignment, main_pdb_array, main_structure):
        lengths_1, lengths_2 = OrderedDict(), OrderedDict()
        missing_from_alignment, missing_from_pdb = [], []
        for i in alignment.template_dict:
            lengths_1[i] = len(alignment.template_dict[i])
            for j in alignment.template_dict[i]:
                try:
                    main_pdb_array[i][j.replace('x','.')]
                except:
                    missing_from_pdb.append(j)
        for i in main_pdb_array:
            lengths_2[i] = len(main_pdb_array[i])
            for j in main_pdb_array[i]:
                try:
                    alignment.template_dict[i][j.replace('.','x')]
                except:
                    try:
                        alignment.template_dict[i][j]
                    except:
                        missing_from_alignment.append(j)
        print('Main structure: ', main_structure)
        print('Alignment dict: ', lengths_1)
        print('Main pdb dict:  ', lengths_2)
        print('Missing from Alignment dict:')
        print(missing_from_alignment)
        print('Missing from Main pdb dict: ')
        print(missing_from_pdb)

    def pdb_pir_mismatch(self, pdb_path, pir_path):
        pdb_struct = PDB.PDBParser(QUIET=True).get_structure('model', pdb_path)[0]
        with open(pir_path, 'r') as pir_file:
            pir_temp_seq = pir_file.readlines()[3]
        pir_temp_seq = pir_temp_seq.replace('/','')
        pdb_dict = OrderedDict()
        print('PDB-PIR mismatches:')
        for chain in pdb_struct:
            for res in chain:
                if PDB.Polypeptide.three_to_one(res.get_resname())!=pir_temp_seq[res.get_id()[1]-1]:
                    print(res, pir_temp_seq[res.get_id()[1]-1])
            for i, r in enumerate(pir_temp_seq):
                if r=='-':
                    try:
                        print(chain[str(i+1)])
                    except:
                        pass

    def force_add_template_to_table(self, table, main_structure, list_of_templates):
        temp = OrderedDict()
        for i,j in table.items():
            if i==main_structure:
                temp[i] = j
                for k in list_of_templates:
                    st = Structure.objects.get(pdb_code__index=k)
                    if st not in table:
                        temp[st] = j-1
            else:
                temp[i] = j
        table = temp
        return table

    # def pdb_pir_mismatch(self, pdb_array, model_sequence):
    #     count = 0
    #     for i,j in pdb_array.items():
    #         for k,l in j.items():
    #             try:
    #                 print(k, l[0].get_parent().get_resname(), model_sequence[count])
    #             except:
    #                 try:
    #                     print(k, l, model_sequence[count])
    #                 except:
    #                     pass
    #             count+=1

