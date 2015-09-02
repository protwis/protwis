from django.shortcuts import render

from construct.models import *
from protein.models import Protein, ProteinConformation
from mutation.models import Mutation

from datetime import datetime


# Create your views here.
def detail(request, slug):
    pc = ProteinConformation.objects.get(protein__entry_name=slug,
        protein__sequence_type__slug='mod')   
 
   # get mutations
    mut = Mutation.objects.filter()



    # get constructs
    c = Construct.objects.get(protein_conformation=pc)
    e = ConstructExpression.objects.filter(construct=c)
#    aux = AuxProtein.objects.get(construct=c)
    exprsn = ConstructExpressionSystem.objects.get()
    sol = ConstructSolubilization.objects.get(construct=c)
    
    
    pur = ConstructPurification.objects.get(construct=c)
    purstep = PurificationStep.objects.filter()
   
    cl = ChemicalList.objects.filter()
    crystalexp = ConstructCrystallization.objects.get(construct=c)
     

    context = {'c':c,'pc': pc,'mut':mut,'exprsn': exprsn,'cl':cl,'sol':sol,'pur':pur,'purstep':purstep,  'crystalexp': crystalexp}
    return render(request,'construct/construct_detail.html',context)

