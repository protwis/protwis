from django.shortcuts import get_object_or_404, render, redirect
from django.views import generic
from django.http import JsonResponse, HttpResponse
from django.db.models import Q, F, Func, Value, Prefetch
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.urls import reverse
from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.http import HttpResponseNotAllowed
from django import forms


from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinSegment
from residue.models import Residue
from structure.models import Structure, StructureModel, StructureExtraProteins
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from common.selection import Selection
from common.views import AbsBrowseSelection
from ligand.models import Ligand, LigandID

import json
from copy import deepcopy
from collections import OrderedDict
import openpyxl
import os


class LandingPage(TemplateView):
    template_name = 'mapper/data_mapper_landing.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context
    
    def post(self, request, *args, **kwargs):

        #############################################################
        ### This method handles POST requests for form submission ###
        #############################################################
        
        if request.method == 'POST':
            
            #################################
            # Utilize ExcelUploadForm class #
            #################################
             
            form = ExcelUploadForm(request.POST,request.FILES)
            
            ####################
            # If form is valid #
            ####################

            if form.is_valid():

                ####################
                # Get cleaned data #
                ####################
                
                file = form.cleaned_data['file']

                ##########################
                # Check if file is .xlsx #
                ##########################

                if not file.name.endswith('.xlsx'):
                    return render(request, self.template_name, {'upload_status': 'Failed'})
                else:

                    ##################################################
                    # Load excel file (workbook) and get sheet names #
                    ##################################################
                    
                    try:
                        workbook = openpyxl.load_workbook(filename=file,read_only=False)

                        sheet_names = workbook.sheetnames

                        ##################################
                        # For each sheet in the workbook #
                        ##################################

                        for sheet_name in sheet_names:
                            
                            ########################
                            # Initialize worksheet #
                            ########################

                            worksheet = workbook.get_sheet_by_name(sheet_name)
                            
                            ###############################
                            # check if headers is correct #
                            ###############################

                            header_list = [worksheet['A1'].value, worksheet['B1'].value, worksheet['C1'].value]
                            
                            ###################################################
                            # If first sheet is receptor with correct headers #
                            ###################################################

                            if sheet_name == 'Phylogenetic Tree' and header_list == ['Receptor (Uniprot)', '1. Feature (Inner cicle)', '2. Order (Outer cicle 1)']:
                                
                                print("success!")
                                
                                #########################################
                                ## if everything is good in do a query ##
                                ##  Fetch: Protein, family and class   ##
                                #########################################

                                protein_data = Protein.objects.filter(entry_name__endswith='_human').prefetch_related('family__parent__parent', 'family__parent__parent__parent').distinct('entry_name')
                                
                                ######################################
                                # Initialize sets for unique entries #
                                ######################################
                                
                                unique_entry_names = set()
                                unique_protein_families = set()
                                unique_protein_classes = set()

                                ##################################
                                # for each entry add to sets of  # 
                                # proteins, families and classes #
                                ##################################
                                
                                for entry in protein_data:
                                    
                                    ########################
                                    # Initiate the entries #
                                    ######################## 
                                    
                                    protein = str(entry)
                                    protein_family = str(entry.family.parent.parent.name)
                                    protein_class = str(entry.family.parent.parent.parent.name)
                                    
                                    ##########################
                                    # Populate the set-lists #
                                    ##########################
                                    
                                    unique_entry_names.add(protein)
                                    unique_protein_families.add(protein_family)
                                    unique_protein_classes.add(protein_class)

                                #########################
                                # convert sets to lists #
                                #########################
                                
                                list_unique_entry_names = list(unique_entry_names)

                                ###################################
                                ## Retrieve all cell values from ##
                                ## columns A, B, and C as lists  ##
                                ###################################
                                
                                ############################
                                # Initiate lists and dicts #
                                ############################

                                headers = [cell.value for cell in worksheet[1]]
                                data_types = [cell.value for cell in worksheet[2]]

                                Correct_values = {}
                                Incorrect_values = {}
                                missing_data_columns = {}

                                ##########################
                                # Run through excel file #
                                ##########################
                                
                                # Identify completely empty columns and skip them
                                empty_columns = set()
                                for col_idx in range(len(headers)):
                                    if all(cell is None for cell in worksheet.iter_cols(min_col=col_idx + 1, max_col=col_idx + 1, min_row=3, values_only=True)):
                                        empty_columns.add(col_idx)
                                print(empty_columns)
                                # Iterate through rows starting from the third row
                                for index, row in enumerate(worksheet.iter_rows(min_row=3, values_only=True), start=3):
                                    # Check the "Receptor (Uniprot)" column for correct values
                                    if row[0] in list_unique_entry_names:
                                        Correct_values[index] = 'Entry correct'
                                    else:
                                        Incorrect_values[index] = 'Wrong entry'

                                    # Check each column for data points, boolean values, and float values
                                    for col_idx, value in enumerate(row):
                                        if col_idx == 0 or col_idx in empty_columns:
                                            continue  # Skip the "Receptor (Uniprot)" column and completely empty columns

                                        if value is not None:
                                            if col_idx not in missing_data_columns:
                                                missing_data_columns[col_idx] = True  # Column has data points
                                            if data_types[col_idx] == 'Boolean' and str(value).lower() not in ['yes', 'no', '1', '0']:
                                                if col_idx not in Incorrect_values:
                                                    Incorrect_values[col_idx] = 'Non-boolean value found'
                                            elif data_types[col_idx] == 'Float':
                                                try:
                                                    float(value)
                                                except ValueError:
                                                    if col_idx not in Incorrect_values:
                                                        Incorrect_values[col_idx] = 'Non-float value found'
                                        else:
                                            if col_idx not in missing_data_columns:
                                                missing_data_columns[col_idx] = False  # Column has missing data points
                                            elif missing_data_columns[col_idx] is True:
                                                missing_data_columns[col_idx] = False  # Column has mixed data points

                                # Prepare the final output
                                for col_idx, has_data in missing_data_columns.items():
                                    if not has_data:
                                        if col_idx not in Incorrect_values:
                                            Incorrect_values[col_idx] = 'Missing data points found'

                                print("Correct Values: ", Correct_values)
                                print("Incorrect Values: ", Incorrect_values)

                                if Incorrect_values:
                                    print("Incorrect values found")
                                    return render(request, self.template_name,{'upload_status': 'Success', 'report_status': 'Failed', 'incorrect_values': Incorrect_values})
                                else:
                                    print("All values are correct")
                                    return render(request, self.template_name,{'upload_status': 'Success', 'report_status': 'Success'})
                    except:
                        return render(request, self.template_name, {'upload_status': 'Failed'})
            else:    
                # Return a 405 Method Not Allowed response if not a POST request
                #return HttpResponseNotAllowed(['POST'])
                return render(request, self.template_name, {'upload_status': 'Failed'})
    
# def LandingPage(request):
#     return render(request, 'mapper/data_mapper_landing.html')

#######################
## Excel upload form ##
#######################

class ExcelUploadForm(forms.Form):
    file = forms.FileField()