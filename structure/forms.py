

from django import forms
import datetime
from django.forms.extras.widgets import SelectDateWidget
from django.forms import ModelForm, Form
from django.utils.safestring import mark_safe

from functools import partial
DateInput = partial(forms.DateInput, {'class': 'datepicker'})
#controlled vocabularies

WT_AA = ( 
    ('H'),
    ('R'),
    ('K'),
    ('D'),
    ('E'),
    ('C'),
    ('G'),
    ('Q'),
    ('N'),
    ('S'),
    ('Y'),
    ('T'),
    ('F'),
    ('A'),
    ('L'),
    ('M'),
    ('I'),
    ('V'),
    ('P'),
    ('W'),
)  


WT_MUT_AA = ( 
    ('H'),
    ('R'),
    ('K'),
    ('D'),
    ('E'),
    ('C'),
    ('G'),
    ('Q'),
    ('N'),
    ('S'),
    ('Y'),
    ('T'),
    ('F'),
    ('A'),
    ('L'),
    ('M'),
    ('I'),
    ('V'),
    ('P'),
    ('W'),
) 

CONSTRUCT_POSITION=(
    ('Please Select'),
    ('N-term'),
    ('Within Receptor'),
    ('C-term'),
    
)  

FUSION_PROTEIN=(
    ('Please Select'),
    ('T4 Lysozyme (T4L)'),
    ('mT4L'),
    ('dsT4L'),
    ('BRIL (b562RIL)'),
    ('Rubredoxin'),
    ('PGS (Pyrococcus abyssi glycogen synthase)'),
     
)

SIGNAL_PEPTIDE=(
    ('Please Select'),
    ('hemagglutinin'),
    ('Other'),
    
)

TAG=(
    ('Please Select'),
    ('FLAG'),
    ('His(10)'),
    ('His(8)'),
    ('His(6)'),
    ('T7'),
    ('TrxA'),
    ('1D4'),
    ('Rho9'),
    ('GFP'),
    ('MBP (Maltose Binding protein)'),
    ('GST'),
    ('Other'),  
)  

PROTEOLYTIC=(
    ('Please Select'),
    ('Thrombin'),
    ('TEV protease'),
    ('HRV 3C protease (PreScission protease)'),
    ('Carboxypeptidase A'),
    ('V8 protease (endoproteinase GluC from Staph8)'),
    ('Other'),
     
)  

PRESENCE=(
    ('Please Select'),
    ('YES'),
    ('NO'),     
)


EXPRESSION_METHOD=(
    ('Select'),
    ('Native_source'),
    ('Stable_cells'),
    ('Transient_transfection'),
    ('Baculovirus'),
    ('Transgenic_organism'),
    ('Other [In case of E.Coli or Yeast recombinant expression]'),
     
) 

HOST_CELL_TYPE=(
    ('Select'),
    ('Sf9'),
    ('Sf21'),
    ('Hi5'),
    ('HEK293'),
    ('HEK293T'),
    ('HEK293EBNA'),
    ('HEK293S(TetR) GnTI−'),
    ('COS1'),
    ('CHOK1'),
    ('ICL3 3rd'),
    ('R. sphaeroides'),
    ('S. cerevisiae'),
    ('E. coli BL21(DE3)'),
    ('P. pastoris SMD1163'),
    ('other [See next field]'),  
)  

HOST_CELL=(
    ('Select'),
    ('Mammal'),
    ('Insect'),
    ('Yeast'),
    ('Bacteria'),
    ('Cell-free'),
    ('other [See next field]'),
     
) 


DETERGENT_TYPE=(
    ('Select'),
    ('DDM'),
    ('DM'),
    ('LMNG'),
    ('OG'),
    ('other [See next field]'),
) 

CHEMICAL_ENZYMATIC_TREATMENT=(
    ('Select'),
    ('Deglycosylation by PNGase F'),
    ('Deglycosylation by Endo F'),
    ('Deglycosylation by Endo H'),
    ('Alkylation by Iodoacetamide'),
    ('Proteolysis by TEV protease'),
    ('Proteolysis by Thrombin'),
    ('Proteolysis by V8 protease'),
    ('Proteolysis by HRV 3C protease'),
    ('Proteolysis by Trypsin'),
    ('Proteolysis by Chymotrypsin'),
    ('Reductive alkylation by Formaldehyde and dimethylaminoborane'),
    ('None'),
    ('P. pastoris SMD1163'),
    ('Other [See remark]'),
      
)  

CRYSTALLIZATION_TYPE=(
    ('Please Select'),
    ('lipidic cubic phase'),
    ('sponge phase'),
    ('bicelles'),
    ('detergent'),
    ('amphipol'),
    ('other [See next field]'),
) 

CRYSTALLIZATION_METHOD=(
    ('Please Select'),
    ('microbatch'),
    ('vapour diffusion'),
    ('free interface diffusion'),
    ('other [See next field]'),
) 


LCP_LIPID=(
    ('Please Select'),
    ('monoolein (9.9 MAG)'),
    ('7.7 MAG'),
    ('7.8 MAG'),
    ('7.9 MAG'),
    ('8.6 MAG'),
    ('8.7 MAG'),
    ('8.8 MAG'),
    ('8.9 MAG'),
    ('monopalmitolein (9.7 MAG)'),
    ('monovaccenin(11.7 MAG)'),
    ('monoeicosenoin (11.9 MAG)'),
    ('other [See next field]'),
) 


CRYSTAL_DETERGENT=(
    ('Please Select'),
    ('LMNG'),
    ('DMNG'),
    ('DGNG'),
    ('OGNG'),
    ('CYMAL5NG'),
    ('CYMAL6NG'),
    ('CYMAL7NG'),
    ('DDM'),
    ('UDM'),
    ('DM'),
    ('NM'),
    ('NG'),
    ('OG'),
    ('HTG'),
    ('HEGA-10'),
    ('LDAO'),
    ('DDAO'),
    ('Cymal5'),
    ('Cymal6'),
    ('Cymal7'),
    ('CHAPS'),
    ('CHAPSO'),   
    ('FacadeEM'),
    ('amphipol A8-35'),
    ('other [See next field]'), 
)  


LIPID=(
    ('Please Select'),
    ('CHS'),
    ('POPC'),
    ('POPE'),
    ('POPS'),
    ('PIP'),
    ('DOPC'),
    ('DMPC'),
    ('DHPC'),
    ('Cholesterol'),
    ('Brain_lipid_extract'),
    ('other [See next field]'),
) 


PROTEIN_TYPE=(
    ('type','Please select Type'),
    ('fusion','Fusion Protein'),
    ('signal','Signal Peptide'),
    ('tag','Tag'),
    ('linker','Linker sequence'),
    ('prot_cleavage','Proteolytic Cleavage site'),  
)


MOD_POS=(
    ('select','Select'),
    ('single','Single'),
    ('pair','Pair'),
    ('range','Range'),  
)

PH=(
    ('single_ph','Single'),
    ('range_ph','Range'),  
)

DELETION=(
    ('select_del','Select'),
    ('del_single','Single'),
    ('del_range','Range'),  
    
)


LIGAND_ID_TYPE=(
    ('Please Select'),
    ('PubChem CID'),
    ('ChEMBL Compound ID'),
    ('SMILES'), 
    ('FASTA sequence (peptide)') ,
    ('UniProt Entry Code (peptide)'),
)


INSERT_POSITION=(
    ('select_ins','Select'),
    ('ins_single','Single'),
    ('ins_range','Range'),  
       
)


class construct_form(forms.Form):

    name_cont = forms.CharField(label='Name', label_suffix='' , max_length=50, required=True)
    #date=forms.DateField(widget = SelectDateWidget) #or alternatively (initial=datetime.date.today)
    pi_name = forms.CharField(label="Name of PI (group leader)", label_suffix='',max_length=50, required=True)
    pi_email=forms.EmailField(label="PI Contact e-mail address", max_length=50,label_suffix='', required=True)  #how to add unique, invalid ?
    url = forms.CharField(initial='http://Database/', widget=forms.TextInput(attrs={'class':'form-control'}) )
    address = forms.CharField(label='Affiliation address', label_suffix='', widget=forms.Textarea, required=True)
    date = forms.DateField(widget=DateInput())
    

#class construct_info(forms.Form):

    pdb = forms.CharField(label='PDB Code', label_suffix='' , max_length=50, required=True)
    pdb_name = forms.CharField(label='Name', label_suffix='' , max_length=50, required=True)
    uniprot=forms.CharField(label_suffix='' , max_length=50, required=True,
    	                    label=mark_safe('Protein name or ID (Uniprot) <a href="http://www.uniprot.org/"</a>(?)'))
    

    ligand_name= forms.CharField(label='Ligand Name', label_suffix='' , max_length=50, required=True)
    ligand_id_type=forms.ChoiceField(choices=[(x,x) for x in LIGAND_ID_TYPE], required=True)
    ligand_activity= forms.CharField(label='Ligand Activity', label_suffix='' , max_length=50, required=True)
    ligand_id= forms.CharField(label='Ligand ID', label_suffix='' , max_length=50, required=True)
    ligand_conc= forms.CharField(label='Ligand Concentration', widget=forms.TextInput(attrs={'class':'half','placeholder': 'value'}), label_suffix='' , max_length=50, required=True)

    deletion_single=forms.CharField(label='Single Deletion', widget=forms.TextInput(attrs={'class':'numeric form-control hidetext del_single del_type row_id', 'placeholder': 'position'}), required=True)
    delet_start = forms.CharField(label='Deletion Start Position',widget=forms.TextInput(attrs={'class':'form-control delet_range hidetext del_range del_type row_id', 'placeholder': 'start'}), required=True)
    delet_end = forms.CharField(label='Deletion End Position', widget=forms.TextInput(attrs={'class':'form-control delet_range hidetext del_range del_type row_id', 'placeholder': 'end'}), required=True)
    deletion=forms.ChoiceField(label='Deletion', choices=DELETION,widget=forms.Select(attrs={'class':'form-control deletion_type half row_id'}), label_suffix='' , required=True)
    
    mutation_details = forms.CharField(label='Mutation Details', label_suffix='' , max_length=50, required=True)
    aa_no = forms.CharField(label='AA No.', widget=forms.TextInput(attrs={'class':'form-control second_table row_id'}), max_length=50, required=True)
    wt_aa=forms.ChoiceField(choices= [(x,x) for x in sorted(WT_AA)], widget=forms.Select(attrs={'class':'form-control mut row_id'}), label='WT AA', label_suffix='' , required=True)
    mut_aa=forms.ChoiceField(choices=[(x,x) for x in sorted(WT_MUT_AA)],widget=forms.Select(attrs={'class':'form-control mut row_id'}),required=True)
    
#class auxiliary_proteins(forms.Form):

    position=forms.ChoiceField(choices=[(x,x) for x in CONSTRUCT_POSITION] , widget=forms.Select(attrs={'class':'custom_select form-control position col_id'}), label='Position in Construct', label_suffix='' , required=True)
    protein_type=forms.ChoiceField(choices=PROTEIN_TYPE, widget=forms.Select(attrs={'class':'custom_select form-control protein_type col_id'}),  label_suffix='' , required=True)
    fusion_prot=forms.ChoiceField(choices=[(x,x) for x in FUSION_PROTEIN], widget=forms.Select(attrs={'class':'custom_select form-control prot_type fusion hidetext col_id'}), label='Fusion Protein', label_suffix='' , required=True)
    
    signal=forms.ChoiceField(choices=[(x,x) for x in SIGNAL_PEPTIDE], widget=forms.Select(attrs={'class':'custom_select form-control prot_type signal hidetext col_id'}), label='signal Peptide', label_suffix='' , required=True)
    other_signal=forms.CharField( widget=forms.TextInput(attrs={'class':'form-control hidetext prot_type others other_signal col_id', 'placeholder':'insert signal'}),required=True)

    tag=forms.ChoiceField(choices=[(x,x) for x in TAG], widget=forms.Select(attrs={'class':'custom_select form-control prot_type tag hidetext col_id'}),label='Tag', label_suffix='' , required=True)
    other_tag=forms.CharField( widget=forms.TextInput(attrs={'class':'hidetext prot_type others other_tag col_id', 'placeholder':'insert tag'}),required=True)


    linker_seq = forms.CharField(label='Linker Sequence', widget=forms.TextInput(attrs={'class':'form-control prot_type linker hidetext col_id', 'placeholder': 'Please input sequence'}), label_suffix='' , max_length=50, required=True)
    prot_cleavage=forms.ChoiceField(choices=[(x,x) for x in PROTEOLYTIC], widget=forms.Select(attrs={'class':'custom_select form-control prot_type prot_cleavage hidetext col_id'}), label='Proteolytic Cleavage Site', label_suffix='' , required=True)
    other_prot_cleavage=forms.CharField( widget=forms.TextInput(attrs={'class':'hidetext others prot_type other_prot col_id', 'placeholder':'insert cleavage site'}),required=True)
    
    insert_pos_type= forms.ChoiceField(label='Insert type', choices=INSERT_POSITION,widget=forms.Select(attrs={'class':'form-control with_rec hidetext half insert_pos_type col_id'}), required=True)
    insert_pos_single=forms.CharField(widget=forms.TextInput(attrs={'class':'form-control with_rec_val hidetext ins_pos_type to_change ins_single custom_single col_id', 'placeholder':'pos'}), label_suffix='' ,  required=True)
    ins_start=forms.CharField(widget=forms.TextInput(attrs={'class':'form-control with_rec_val hidetext ins_start ins_pos_type ins_range custom_range col_id', 'placeholder':'pos'}), label_suffix='' ,  required=True)
    ins_end=forms.CharField(widget=forms.TextInput(attrs={'class':'form-control with_rec_val hidetext ins_end ins_pos_type ins_range custom_range col_id', 'placeholder':'pos'}), label_suffix='' ,  required=True)
    

    presence=forms.ChoiceField(choices=[(x,x) for x in PRESENCE], widget=forms.Select(attrs={'class':'custom_select form-control col_id col_id'}) ,required=True)
    

#class modification(forms.Form):

    aamod=forms.CharField(label='aamod', widget=forms.TextInput(attrs={'class':'form-control widen second_table row_id' }), max_length=50, required=True)
    aamod_position=forms.ChoiceField(choices=MOD_POS,widget=forms.Select(attrs={'class':'form-control half aamod_pos_type row_id'}),label='aamod position', label_suffix='' , required=True)
    
    aamod_single=forms.CharField(label='aamod_pos',widget=forms.TextInput(attrs={'class':'form-control hidetext aa_type custom_single single row_id', 'placeholder':'insert position'}), label_suffix='' , max_length=50, required=True)
    aamod_pair_one=forms.CharField(label='aamod_pos',widget=forms.TextInput(attrs={'class':'form-control hidetext aa_type custom_pair pair row_id', 'placeholder':'pos'}), label_suffix='' , max_length=50, required=True)
    aamod_pair_two=forms.CharField(label='aamod_pos',widget=forms.TextInput(attrs={'class':'form-control hidetext aa_type custom_pair pair row_id', 'placeholder':'pos'}), label_suffix='' , max_length=50, required=True)
    aamod_start=forms.CharField(label='aamod_range',widget=forms.TextInput(attrs={'class':'form-control hidetext aa_type custom_range range row_id', 'placeholder':'pos'}), label_suffix='' , max_length=50, required=True)
    aamod_end=forms.CharField(label='aamod_range',widget=forms.TextInput(attrs={'class':'form-control hidetext aa_type custom_range range row_id', 'placeholder':'pos'}), label_suffix='' , max_length=50, required=True)
    mod_remark = forms.CharField(widget=forms.Textarea(attrs={'class':'form-control modremark modremark row_id'}), required=True)

#class expression(forms.Form):

    expr_method=forms.ChoiceField(choices=[(x,x) for x in EXPRESSION_METHOD], label='Expression Method')
    host_cell_type=forms.ChoiceField(choices=[(x,x) for x in HOST_CELL_TYPE], label='Host Cell Type', label_suffix='' , required=True)
    host_cell=forms.ChoiceField(choices=[(x,x) for x in HOST_CELL], label='Host Cell', label_suffix='' , required=True)
    expr_remark = forms.CharField(label='Remark', label_suffix='', widget=forms.Textarea,required=True)


#class solubilization(forms.Form):
    deterg_type=forms.ChoiceField(choices=[(x,x) for x in DETERGENT_TYPE], label='Detergent Type', label_suffix='' , required=True)
    other_deterg_type=forms.CharField( widget=forms.TextInput(attrs={'class':'hidetext other_type_deterg', 'placeholder':'insert type'}),label='other deterg type', max_length=50, required=True)

    deterg_concentr = forms.CharField(label='Detergent Concentration', label_suffix='' ,widget=forms.TextInput(attrs={'class':'form-control half', 'placeholder': 'value'}), max_length=50, required=True) 
    deterg_concentr_unit = forms.CharField(label='Detergent Concentration Unit', widget=forms.TextInput(attrs={'class':'form-control unit','placeholder': 'unit i.e.% w/v'}), label_suffix='' , max_length=50, required=True)  
    solub_additive = forms.CharField(label='Solubilization additive', label_suffix='' , max_length=50, required=True)
    additive_concentr= forms.CharField(label='Additive concentration', widget=forms.TextInput(attrs={'class':'form-control half', 'placeholder':'value'}),label_suffix='' , max_length=50, required=True)
    addit_concentr_unit = forms.CharField(label='Additive Concetration Unit', widget=forms.TextInput(attrs={'class':'form-control unit','placeholder': 'unit i.e.:% w/v'}), label_suffix='' , max_length=50, required=True)  
    
    chem_enz_treatment=forms.ChoiceField(choices=[(x,x) for x in CHEMICAL_ENZYMATIC_TREATMENT],widget=forms.Select(attrs={'class':'form-control chem row_id'}) , label='Chemical Enzymatic Treatment', label_suffix='' , required=True)
    sol_remark = forms.CharField(label='Remark', label_suffix='', widget=forms.TextInput(attrs={'class':'form-control hidetext chem_enz_remark row_id', 'placeholder':'insert treatment'}),required=True)

#class crystal_conditions(forms.Form):
    crystal_type=forms.ChoiceField(choices=[(x,x) for x in CRYSTALLIZATION_TYPE], label='Type', label_suffix='' , required=True)
    other_crystal_type= forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext other_cryst_type', 'placeholder':'insert type'}))    

    crystal_method=forms.ChoiceField(choices=[(x,x) for x in CRYSTALLIZATION_METHOD], label='Method', label_suffix='' , required=True)
    other_method=forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext other_method', 'placeholder':'insert method'}),  max_length=50, required=True)
    protein_concentr= forms.CharField(label='Protein Concentration',widget=forms.TextInput(attrs={'class':'half', 'placeholder':'value'}), label_suffix='' , max_length=50, required=True)
    protein_conc_unit= forms.CharField(label='Protein Concentration Unit', widget=forms.TextInput(attrs={'class':'unit','placeholder': 'unit i.e.mg/ml'}), label_suffix='' , max_length=50, required=True)
    temperature= forms.CharField(label='Temperature', label_suffix='' , max_length=50,widget=forms.TextInput(attrs={'placeholder': '°C', 'type':'numberzzz'}), required=True)
    
    ph= forms.ChoiceField(choices=PH,widget=forms.Select(attrs={'class':'ph half'}), required=True)
    ph_single=forms.CharField(widget=forms.TextInput(attrs={'class':'ph_type custom_single single_ph', 'placeholder':'insert pH'}))
    ph_range_one=forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext ph_type custom_range range_ph', 'placeholder':'pH'}))
    ph_range_two=forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext ph_type custom_range range_ph', 'placeholder':'pH'}))
   
    crystal_remark = forms.CharField(label='Remark', label_suffix='', widget=forms.Textarea,required=True)
    lcp_lipid=forms.ChoiceField(choices=[(x,x) for x in LCP_LIPID], label='LCP Lipid', label_suffix='' , required=True)
    other_lcp_lipid= forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext other_lcp', 'placeholder':'insert lipid'}))

    lcp_add= forms.CharField(label='LCP Lipid Additive', label_suffix='' , max_length=50, required=True)
    lcp_conc= forms.CharField(label='LCP Lipid Additive Concentration', widget=forms.TextInput(attrs={'class':'half','placeholder': 'value'}), label_suffix='' , max_length=50, required=True)
    lcp_conc_unit=forms.CharField(label='LCP Lipid Additive Concentration', widget=forms.TextInput(attrs={'class':'unit','placeholder': 'unit'}), label_suffix='' , max_length=50, required=True)
    detergent=forms.ChoiceField(choices=[(x,x) for x in CRYSTAL_DETERGENT], label='Detergent', label_suffix='' , required=True)
    deterg_conc= forms.CharField(label='Deterg Concentration', widget=forms.TextInput(attrs={'class':'half','placeholder': 'value'}), label_suffix='' , max_length=50, required=True)
    
    deterg_conc_unit=forms.CharField(label='Det conc unit', widget=forms.TextInput(attrs={'class':'unit','placeholder': 'unit i.e.%w/v'}), label_suffix='' , max_length=50, required=True)
    other_deterg=forms.CharField(label='other deterg', widget=forms.TextInput(attrs={'class':'hidetext other_det', 'placeholder':'insert detergent'}), label_suffix='' , max_length=50, required=True)

    lipid=forms.ChoiceField(choices=[(x,x) for x in LIPID], label='Lipid', label_suffix='' , required=True)
    other_lipid=forms.CharField(widget=forms.TextInput(attrs={'class':'hidetext other_lipid', 'placeholder':'insert lipid'}),  max_length=50, required=True)

    lipid_concentr= forms.CharField(label='Lipid Concentration', widget=forms.TextInput(attrs={'class':'half', 'placeholder':'value'}),label_suffix='' , max_length=50, required=True)
    lipid_concentr_unit=forms.CharField(label='Lipid conc unit', widget=forms.TextInput(attrs={'class':'unit','placeholder': 'unit i.e.%w/v'}), label_suffix='' , max_length=50, required=True)

    ligand_conc_unit=forms.CharField(widget=forms.TextInput(attrs={'class':'unit','placeholder': 'unit i.e.%w/v'}))
    
    chemical_comp= forms.CharField(label='Chemical Component',widget=forms.TextInput(attrs={'class':'form-control row_id row_id'}),  max_length=50, required=True)
    concentr= forms.CharField(label='Concentration',  widget=forms.TextInput(attrs={'class':'form-control half half row_id', 'placeholder':'value'}), label_suffix='' , max_length=50, required=True)
    concentr_unit=forms.CharField(label='Concentration unit', widget=forms.TextInput(attrs={'class':'form-control unit unit row_id','placeholder': 'unit'}), label_suffix='' , max_length=50, required=True)
    
    #concentr_unit.widget.attrs.update({'class' : 'form-control'})

f = construct_form()


