from django.shortcuts import render

from construct.models import *
from protein.models import Protein, ProteinConformation
from mutation.models import Mutation

from datetime import datetime


# Create your views here.
def detail(request, slug):
#    construct=Construct.objects.all()
    pc = ProteinConformation.objects.get(protein__entry_name=slug,
        protein__sequence_type__slug='mod')   
 
   # get mutations
    mut = Mutation.objects.filter()



    # get constructs
    c = Construct.objects.get(protein_conformation=pc)
    e = ConstructExpression.objects.filter(construct=c)
    exprsn = ConstructExpressionSystem.objects.get()
    sol = ConstructSolubilization.objects.get(construct=c)
    
    
    pur = ConstructPurification.objects.get(construct=c)
    purstep = PurificationStep.objects.filter()
#    pursteptype = PurificationStepType.objects.get(purificationstep=purstep)
   
    cl = ChemicalList.objects.filter()
    crystalexp = ConstructCrystallization.objects.get(construct=c)
     

#    context = {'p': p}
    context = {'c':c,'pc': pc,'mut':mut, 'exprsn': exprsn,'cl':cl,'sol':sol,'pur':pur,'purstep':purstep,  'crystalexp': crystalexp}
#    print (crystalexp)
    return render(request,'construct/construct_detail.html',context)

