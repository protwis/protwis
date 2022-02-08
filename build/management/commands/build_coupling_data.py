import csv
import logging
import os
import sys
from itertools import islice
from collections import defaultdict

import xlrd

from common.models import Publication, WebLink, WebResource
from django.conf import settings
from django.core.management.base import BaseCommand

from protein.models import Protein, ProteinFamily, ProteinCouplings
from ligand.models import Ligand
from ligand.functions import get_or_make_ligand

class Command(BaseCommand):
    help = 'Build G proteins'

    # source files
    iupharcoupling_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'iuphar_coupling_data.csv'])
    master_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'GPCR-G_protein_couplings.xlsx'])

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        self.purge_coupling_data()
        self.logger.info('PASS: purge_coupling_data')
        self.create_iuphar_couplings()
        self.logger.info('PASS: create_iuphar_couplings')
        self.create_data_couplings()
        self.logger.info('PASS: create_data_couplings')

    def purge_coupling_data(self):
        """DROP data from the protein_gprotein_pair table."""
        try:
            ProteinCouplings.objects.filter().delete()
        except Exception as msg:
            self.logger.warning('Existing protein couplings cannot be deleted' + str(msg))

    def create_iuphar_couplings(self, filenames=False):
        """
        Function to add IUPHAR receptor-G protein couplings.

        The function reads a iupharcoupling_file, which comes from parsing the Guide_to_Pharmacology.
        """
        self.logger.info('CREATING: IUPHAR couplings')
        print("PROCESSING: IUPHAR couplings")

        translation = {'Gs family': ['Gs', '100_001_001'],
                       'Gi/Go family': ['Gi/o', '100_001_002'],
                       'Gq/G11 family': ['Gq/11', '100_001_003'],
                       'G12/G13 family': ['G12/13', '100_001_004'],
                       'GPa1 family': ['GPa1 family', '100_001_005']
                       }

        filenames = [self.iupharcoupling_file]

        source = "GuideToPharma"
        for filename in filenames:
            filepath = self.iupharcoupling_file
            self.logger.info('Reading filename ' + filename)
            pub_years = defaultdict(int)
            pub_years_protein = defaultdict(set)

            with open(filepath, 'r') as f:
                reader = csv.reader(f)
                for row in islice(reader, 1, None):  # skip first line
                    entry_name = row[3]
                    primary = row[11]
                    secondary = row[12]
                    primary_pubmed = row[15]
                    secondary_pubmed = row[16]
                    # fetch protein
                    try:
                        p = Protein.objects.get(entry_name=entry_name)
                    except Protein.DoesNotExist:
                        self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                        print("protein not found for ", entry_name)
                        continue

                    primary = primary.replace("G protein (identity unknown)", "None")  # replace none
                    primary = primary.split(", ")

                    secondary = secondary.replace("G protein (identity unknown)", "None")  # replace none
                    secondary = secondary.split(", ")

                    if primary == 'None' and secondary == 'None':
                        print('no data for ', entry_name)
                        continue

                    try:
                        for gp in primary:
                            if gp in ['', 'None', '_-arrestin', 'Arrestin',
                                      'G protein independent mechanism']:  # skip bad ones
                                continue
                            g = ProteinFamily.objects.get_or_create(name=translation[gp][0], slug=translation[gp][1])[0]
                            # print(p, g)
                            gpair = ProteinCouplings(protein=p, g_protein=g, transduction='primary', source=source)
                            gpair.save()

                            for pmid in primary_pubmed.split("|"):
                                try:
                                    test = int(pmid)
                                except:
                                    continue
                                try:
                                    pub = Publication.objects.get(web_link__index=pmid,
                                                                  web_link__web_resource__slug='pubmed')
                                except Publication.DoesNotExist as e:
                                    pub = Publication()
                                    try:
                                        pub.web_link = WebLink.objects.get(index=pmid,
                                                                           web_resource__slug='pubmed')
                                    except WebLink.DoesNotExist:
                                        wl = WebLink.objects.create(index=pmid,
                                                                    web_resource=WebResource.objects.get(slug='pubmed'))
                                        pub.web_link = wl
                                pub.update_from_pubmed_data(index=pmid)
                                pub.save()
                                pub_years[pub.year] += 1
                                pub_years_protein[pub.year].add(entry_name)
                                gpair.references.add(pub)

                    except Exception as e:
                        print("error in primary assignment", p, gp, e)

                    try:
                        for gp in secondary:
                            if gp in ['None', '_-arrestin', 'Arrestin', 'G protein independent mechanism',
                                      '']:  # skip bad ones
                                continue
                            if gp in primary:  # sip those that were already primary
                                continue
                            g = ProteinFamily.objects.get_or_create(name=translation[gp][0], slug=translation[gp][1])[0]
                            gpair = ProteinCouplings(protein=p, g_protein=g, transduction='secondary', source=source)
                            gpair.save()

                            for pmid in secondary_pubmed.split("|"):
                                try:
                                    test = int(pmid)
                                except:
                                    continue

                                try:
                                    pub = Publication.objects.get(web_link__index=pmid,
                                                                  web_link__web_resource__slug='pubmed')
                                except Publication.DoesNotExist as e:
                                    pub = Publication()
                                    try:
                                        pub.web_link = WebLink.objects.get(index=pmid,
                                                                           web_resource__slug='pubmed')
                                    except WebLink.DoesNotExist:
                                        wl = WebLink.objects.create(index=pmid,
                                                                    web_resource=WebResource.objects.get(slug='pubmed'))
                                        pub.web_link = wl
                                pub.update_from_pubmed_data(index=pmid)
                                pub.save()
                                pub_years[pub.year] += 1
                                pub_years_protein[pub.year].add(entry_name)
                                gpair.references.add(pub)
                    except Exception as e:
                        print("error in secondary assignment", p, gp, e)


        self.logger.info('COMPLETED CREATING IUPHAR couplings')

    def assess_ligand_id(self, spreadsheet):
        """
        Function that fetches the ligands data from the master file for the different sources
        and compare them to the LigandIDs sheet. A dictionary is created that contains the information
        to be further parsed via the get_or_make_ligand function. Data from Bouvier and Inoue, processed
        by David Gloriam.
        """
        #rading data by the ligand sheets
        Bligands = spreadsheet.sheet_by_name("B-Ligands")
        Iligands = spreadsheet.sheet_by_name("I-Ligands")
        IDs = spreadsheet.sheet_by_name("LigandIDs")
        #counting sheets rows
        Blig_rows = Bligands.nrows
        Ilig_rows = Iligands.nrows
        IDs_rows = IDs.nrows

        #third value is the column number bearing the ligand name
        #bouvier --> ligand renamed - Bouvier (9)
        #inoue   --> ligand short (author)    (2)
        tupled_sources = [(Bligands, Blig_rows, 9),(Iligands, Ilig_rows, 2)]

        def parseValue(s, type='int'):
            """
            Function to return an int or NA from the ligand data values
            """
            if s == '':
                return 'NA'
            else:
                if type == 'int':
                    return int(s)
                else:
                    return str(s)

        IDsdata = {}
        # sources_data = {}
        #create a dictionary of the complete info in the LigandIDs
        #sheet, using GtP - Ligand ID column as key and setting the
        #GtP - PubChem CID, GtP - Ligand Name, GtP - SMILES, Accession Number
        #total 5 values
        for i in range(1, IDs_rows): #skip the header row
            if IDs.cell_value(i, 4) not in IDsdata.keys():
                IDsdata[IDs.cell_value(i, 4)] = [IDs.cell_value(i, 3), parseValue(IDs.cell_value(i, 7)), parseValue(IDs.cell_value(i, 11), 'str'), IDs.cell_value(i, 1)]

        #Now need to parse the Bouvier and Inoue ligand data
        #search for matches and provide expanded info for get_or_make_ligand function
        for triple in tupled_sources:
            for i in range(1, triple[1]):
                if triple[0].cell_value(i, 5) not in IDsdata.keys():
                    #if the ligand is not in the IDs sheet, save the UniProt and the ligand name
                    IDsdata[triple[0].cell_value(i, 5)] = [triple[0].cell_value(i, 0), triple[0].cell_value(i, triple[2])]

        #this creates a dictionary with all the ligands data from the sources
        #with PubChem CID, SMILES and Ligand Name where available.
        #Otherwise we have the ligand name from the sheet
        return IDsdata


    def assess_variants(self, dictionary, iterator, column_variants):
        """
        Function that assesses the variants of the G-Protein coupling or Arrestin coupling.
        For now variants assessed are: Isoform 2 of Go and Arrestin without GRK, but more can be further
        implemented, thus we need to keep this function highly flexible.
        """
        if iterator in column_variants:
            dictionary['variant'] = column_variants[iterator]
        else:
            dictionary['variant'] = 'Regular'

    def assess_type(self, accession_id):
        """
        Function that assess the type of the ligand tested making a call to the Protein model.
        It is used to provide precise info to the get_or_make_ligand function in terms of
        ligand_type for the loaded properties (ligand_properities table). Usually used for Peptide
        or Protein calls. Defaults as Peptide call.
        """
        call = list(Protein.objects.filter(accession=accession_id).values_list("family__parent__parent__name"))
        label = call[0][0].split(' ')[0].lower()
        if (label != 'peptide') and (label != 'protein'):
            label = 'peptide'
        return label

# this is the new function that will overwrite the function read_coupling
# reading all different sheets of the single file
    def read_all_coupling(self, filenames=False):
        """
        Yet another function to read G-protein coupling data coming in Excel files.
        The idea is that now the format will hopefully be fixed in the same way for data
        coming from different groups. For now the data comes from Bouvier, Inoue and Roth
        and has been processed by David Gloriam.
        """
        book = xlrd.open_workbook(filenames)

        bouvier = book.sheet_by_name("B-import")
        inoue = book.sheet_by_name("I-import")
        roth = book.sheet_by_name("R-import")
        B_rows = bouvier.nrows
        I_rows = inoue.nrows
        R_rows = roth.nrows

        combinations = [(bouvier, B_rows, 'Bouvier'), (inoue, I_rows, 'Inoue'), (roth, R_rows, 'Roth')]
        #variable that sets the number of subunits parsed
        #in the spreadsheet. If more subunits are introduced
        #please update this number
        subunits = 17
        #each key if the section of data in the spreadsheet
        #values are: start column, ArrB2_NO_GRK column, GoA column and GoB column
        variants = ["isoform 1 (GoA)", "isoform 2 (GoB)", "no GRK"]
        variant_indices = [31, 32, 40]
        start_index = 26
        columns = {
                    'logmaxec50':   {"start": start_index, "variants": dict(zip(variant_indices, variants))},
                    'pec50deg':     {"start": start_index+subunits*1, "variants": dict(zip([i+subunits*1 for i in variant_indices], variants))},
                    'emaxdeg':      {"start": start_index+subunits*2, "variants": dict(zip([i+subunits*2 for i in variant_indices], variants))},
                    'stddeg':       {"start": start_index+subunits*3, "variants": dict(zip([i+subunits*3 for i in variant_indices], variants))}
                    }

        data = {}
        """data is a dictionary and must have a format:
        {'<protein>':
            {'<gproteinsubunit>':
                {'logmaxec50': <logmaxec50>,
                 'pec50deg': <pec50deg>,
                 'emaxdeg': <emaxdeg>,
                 'stddeg': <stddeg>}
            }
        }"""

        ligands = self.assess_ligand_id(book)

        #for each sheet of data from souces
        for tuple in combinations:
            data[tuple[2]] = {}
            #for each row, set protein as dict key
            #and fetch fixed columns info
            for i in range(2, tuple[1]):
                protein = tuple[0].cell_value(i, 0)
                protein_dict = {}
                protein_dict['ligand_id'] = tuple[0].cell_value(i, 4)
                if tuple[0].cell_value(i, 5) == 'Physiological':
                    protein_dict['ligand_physiological'] = True
                else:
                    protein_dict['ligand_physiological'] = False
                #for each block of data parse the sheet
                #and retrieve the associated info
                for key in columns.keys():
                    for j in range(columns[key]["start"], (columns[key]["start"] + subunits)):
                        gproteinsubunit = tuple[0].cell_value(1, j).split("\n")[-1]
                        if gproteinsubunit not in protein_dict.keys():
                            protein_dict[gproteinsubunit] = {}
                        textType = tuple[0].cell(i, j).ctype
                        if textType == 5:
                            protein_dict[gproteinsubunit][key] = None
                        else:
                            protein_dict[gproteinsubunit][key] = tuple[0].cell_value(i, j)
                        self.assess_variants(protein_dict[gproteinsubunit], j, columns[key]["variants"])
                #apply temporary dict to master dict
                data[tuple[2]][protein] = protein_dict

        return data, ligands

    def create_data_couplings(self):
        """This function adds all coupling data coming from the master Excel file."""
        self.logger.info('CREATE data couplings')

        # read source files
        filepath = self.master_file
        self.logger.info('Reading file ' + filepath)
        #new function for reading data from the master spreadsheet
        data, ligands = self.read_all_coupling(filepath)

        lookup = {}
        bulk = []
        ######### MODIFIED #########
        for source in data.keys():
            print("PROCESSING: "+str(source)+" DATA")
            for entry_name, couplings in data[source].items():
                # if it has / then pick first, since it gets same protein
                entry_name = entry_name.split("/")[0]
                # print("PROCESSING: "+str(entry_name)+" ROW")
                # Fetch protein
                try:
                    p = Protein.objects.get(entry_name=entry_name.lower()+"_human")
                except Protein.DoesNotExist:
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    print("protein not found for ", entry_name)
                    continue

                gproteins = list(couplings.keys())[2:]

                for header in gproteins:
                    # print("PROCESSING: "+str(header)+" COLUMN")
                    gprotein = header.split('/')[0]
                    if gprotein not in lookup:
                        gp = Protein.objects.filter(family__name=gprotein, species__common_name="Human")[0]
                        lookup[gprotein] = gp
                    else:
                        gp = lookup[gprotein]

                    # Assume they are there.
                    if gp.family.slug not in lookup:
                        # print("SEARCHING FOR FAMILY OF SLUG: " + str("_".join(gp.family.slug.split("_")[:3])))
                        lookup[gp.family.slug] = ProteinFamily.objects.get(slug="_".join(gp.family.slug.split("_")[:3]))

                    g = lookup[gp.family.slug]
                    # print("PROCESSING: LIGANDS INFORMATION OF " + str(ligands[couplings['ligand_id']][0]) + " FOR "+ str(g))
                    #ligand here should be fetched via the function get_or_make_ligand
                    #(0) Ligand Name, (1) PubChem CID, (2) SMILES, (3) Accession Number
                    if len(ligands[couplings['ligand_id']]) == 4:
                        if ligands[couplings['ligand_id']][2] != 'NA':
                            # print("FETCHING: LIGAND " + str(ligands[couplings['ligand_id']][0]) + " BY SMILES")
                            try:
                                l = get_or_make_ligand(ligands[couplings['ligand_id']][2], 'SMILES', ligands[couplings['ligand_id']][0])
                            except UnboundLocalError:
                                # print("ERROR WITH SMILES. TRYING WITH CID")
                                # print("FETCHING: LIGAND " + str(ligands[couplings['ligand_id']][0]) + " BY CID")
                                l = get_or_make_ligand(ligands[couplings['ligand_id']][1], 'PubChem CID', ligands[couplings['ligand_id']][0])
                        elif ligands[couplings['ligand_id']][1] != 'NA':
                            # print("FETCHING: LIGAND " + str(ligands[couplings['ligand_id']][0]) + " BY CID")
                            l = get_or_make_ligand(ligands[couplings['ligand_id']][1], 'PubChem CID', ligands[couplings['ligand_id']][0])
                        else:
                            #make the call to search for peptide/protein
                            # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                            label = self.assess_type(ligands[couplings['ligand_id']][3])
                            # print("ADDING: LIGAND " + str(ligands[couplings['ligand_id']][0]) + " AS " + str(label))
                            l = get_or_make_ligand('NA', 'NA', ligands[couplings['ligand_id']][0], label)
                    else:
                        #make the call to search for peptide/protein
                        # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                        label = self.assess_type(ligands[couplings['ligand_id']][1])
                        # print("ADDING: LIGAND " + str(ligands[couplings['ligand_id']][0]) + " AS " + str(label))
                        l = get_or_make_ligand('NA', 'NA', ligands[couplings['ligand_id']][0], label)

                    #here it fills the data in ProteinCouplings model
                    #new data in dictionary: {'ligand_name': 'Serotonin', 'ligand_id': 5.0, 'ligand_type': 'Physiological'}
                    # print("CREATING CALL TO MODEL")
                    if (couplings[header]['logmaxec50'] != None) or (couplings[header]['pec50deg'] != None) or (couplings[header]['emaxdeg'] != None) or (couplings[header]['stddeg'] != None):
                        gpair = ProteinCouplings(protein=p,
                                                   g_protein=g,
                                                   ligand=l,
                                                   variant=couplings[header]['variant'],
                                                   source=source,
                                                   logmaxec50=couplings[header]['logmaxec50'],
                                                   pec50=couplings[header]['pec50deg'],
                                                   emax=couplings[header]['emaxdeg'],
                                                   stand_dev=couplings[header]['stddeg'],
                                                   physiological_ligand=couplings['ligand_physiological'],
                                                   g_protein_subunit=gp)
                        # print("APPENDING CALL TO BULK")
                        bulk.append(gpair)

        ProteinCouplings.objects.bulk_create(bulk)

        self.logger.info('COMPLETED CREATE data couplings')
