from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *
from django.conf import settings
from django.db.models import Prefetch
from django.utils.text import slugify

from ligand.models import Ligand, LigandType, AssayExperiment, LigandVendors, LigandVendorLink
from protein.models import Protein

import datamol as dm
import datetime
import pandas as pd

class Command(BaseBuild):
    help = "Import ChEMBL data from data repository"
    bulk_size = 50000

    def add_arguments(self, parser):
        parser.add_argument("--test_run",
                            action="store_true",
                            help="Skip this during a test run",
                            default=False)
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing ligand records')

    def handle(self, *args, **options):
        if options["test_run"]:
            print("Skipping in test run")
            return

        if options['purge']:
            # delete any existing ChEMBL bioactivity data
            # For deleting the ligands - purse all ligand data using the GtP ligand build
            print("Started purging ChEMBL bioactivity data")
            self.purge_chembl_data()
            print("Ended purging ChEMBL bioactivity data")

        # Parse ChEMBL ligand data
        print("Started building ChEMBL ligands")
        self.build_chembl_ligands()
        print("Ended building ChEMBL ligands")

        # Parse ChEMBL bioactivity data
        print("\n\nStarted building ChEMBL bioactivities")
        self.build_chembl_bioactivities()
        print("Ended building ChEMBL bioactivities")

        # Parse ChEMBL/PubChem vendor data
        print("\n\nStarted building PubChem vendor data")
        self.build_pubchem_vendor_links()
        print("Ended building PubChem vendor data")



    @staticmethod
    def purge_chembl_data():
        Ligand.objects.all().delete()

    @staticmethod
    def build_chembl_ligands():
        print("\n===============\n#1 Reading ChEMBL ligand data", datetime.datetime.now())
        ligand_input_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_cpds.csv"])
        ligand_data = pd.read_csv(ligand_input_file, keep_default_na=False)
        for column in ligand_data:
            ligand_data[column] = ligand_data[column].replace({"":None})
        print("Found", len(ligand_data), "ligands")

        # Collect ChEMBL IDs from existing ligands
        print("\n#2 Collecting ChEMBL IDs from existing ligands", datetime.datetime.now())
        wr_chembl = WebResource.objects.get(slug="chembl_ligand")
        wr_pubchem = WebResource.objects.get(slug="pubchem")

        # Removing existing ligands based on ChEMBL IDs => second filter in loop necessary because of concatenated non-parent IDs
        existing_ids = list(LigandID.objects.filter(web_resource=wr_chembl).values_list("index", flat=True).distinct())
        filtered_ligands = ligand_data.loc[~ligand_data["molecule_chembl_id"].isin(existing_ids)].copy()

        # Set ChEMBL ID as name for ligands without pref_name
        filtered_ligands.loc[filtered_ligands['pref_name'].isnull(), 'pref_name'] = filtered_ligands.loc[filtered_ligands['pref_name'].isnull(), 'molecule_chembl_id']

        # Parse all new small-molecule ChEMBL ligands
        print("\n#3 Building new small-molecule ChEMBL ligands", datetime.datetime.now())
        sm_data = filtered_ligands.loc[filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide", None])].reset_index()
        lig_entries = len(sm_data)
        print("Found", lig_entries, "new small molecules")

        # Additional matching via PubChem and InchiKeys
        existing_cids = list(LigandID.objects.filter(web_resource=wr_pubchem).values_list("index", flat=True).distinct())
        existing_inchis = list(Ligand.objects.exclude(inchikey=None).values_list("inchikey", flat=True).distinct())

        smallmol = LigandType.objects.get(slug="small-molecule")
        ligands = []
        weblinks = []
        for index, row in sm_data.iterrows():
            insert = True
            ids = [row['molecule_chembl_id']]

            # Filtering on non-parent ChEMBL IDs
            if row['other_ids'] != None:
                extra_ids = row['other_ids'].split(";")
                existing = list(set.intersection(set(extra_ids), set(existing_ids)))
                if len(existing) > 0:
                    # Add missing link to parent ChEMBL ID
                    match = LigandID.objects.get(index=existing[0], web_resource=wr_chembl)
                    LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand_id = match.ligand_id).save()
                    print("Found existing non-parent ChEMBL", existing[0], "for parent", row['molecule_chembl_id'])
                    insert = False # use this switch in case this is the last one in the list skipping the bulk insert
                else:
                    ids = ids + extra_ids

            # Filtering on PubChem CIDs
            if row['pubchem_cid'] != None:
                cids = row['pubchem_cid'].split(";")
                existing = list(set.intersection(set(cids), set(existing_cids)))
                if len(existing) > 0:
                    # Add missing link to parent ChEMBL ID
                    match = LigandID.objects.get(index=existing[0], web_resource=wr_pubchem)
                    LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand_id = match.ligand_id).save()
                    insert = False # use this switch in case this is the last one in the list skipping the bulk insert

            # For those rare cases in which neither the ChEMBL ID nor the PubChem ID was matched, but the InchiKey was
            if insert and row['standard_inchi_key'] in existing_inchis:
                ligand = Ligand.objects.get(inchikey=row['standard_inchi_key'])
                LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand = ligand).save()
                if row['pubchem_cid'] != None:
                    cids = row['pubchem_cid'].split(";")
                    for cid in cids:
                        LigandID(index=cid, web_resource=wr_pubchem, ligand = ligand).save()
                insert = False

            if insert:
                # creating ligand
                ligands.append(Ligand(name = row['pref_name'], ambiguous_alias = False))
                ligands[-1].ligand_type = smallmol
                ligands[-1].smiles = row['smiles']
                ligands[-1].inchikey = row['standard_inchi_key']
                ligands[-1].sequence = row["sequence"]

                try:
                    input_mol = dm.to_mol(row['smiles'], sanitize=True)
                    if input_mol:
                        # Check if InChIKey has been set
                        if ligands[-1].inchikey == None:
                            ligands[-1].inchikey = dm.to_inchikey(input_mol)

                        # Cleaned InChIKey has been set
                        # ligands[-1].clean_inchikey = get_cleaned_inchikey(row['smiles'])

                        # Calculate RDkit properties
                        ligands[-1].mw = dm.descriptors.mw(input_mol)
                        ligands[-1].rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
                        ligands[-1].hacc = dm.descriptors.n_hba(input_mol)
                        ligands[-1].hdon = dm.descriptors.n_hbd(input_mol)
                        ligands[-1].logp = dm.descriptors.clogp(input_mol)
                except:
                    skip = True

                # Adding ligand IDs
                for id in ids:
                    weblinks.append({"link" : LigandID(index=id, web_resource=wr_chembl), "lig_idx" : len(ligands)-1})
                if row['pubchem_cid'] != None:
                    cids = row['pubchem_cid'].split(";")
                    for cid in cids:
                        weblinks.append({"link" : LigandID(index=cid, web_resource=wr_pubchem), "lig_idx" : len(ligands)-1})

            # BULK insert every X entries or last entry
            if len(ligands) == Command.bulk_size or (index == lig_entries - 1):
                Ligand.objects.bulk_create(ligands)

                # Ligands have been inserted => update LigandIDs for pairing
                for pair in weblinks:
                    pair["link"].ligand = ligands[pair["lig_idx"]]
                LigandID.objects.bulk_create([pair["link"] for pair in weblinks])

                print("Inserted", index + 1, "out of", lig_entries, "ligands")
                ligands = []
                weblinks = []

        # Parse all new non-small-molecule ChEMBL ligands
        print("\n#4 Building new non-small-molecule ChEMBL ligands", datetime.datetime.now())
        nonsm_data = filtered_ligands.loc[~filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide", None])]
        print("Found", len(nonsm_data), "new non-small-molecules")

        ligands = []
        ligand_types = {"Unknown": "na", "Protein": "protein"}
        weblinks = []
        for index, row in nonsm_data.iterrows():
            ids = {}
            ids["smiles"] = row['smiles']
            ids["sequence"] = row['sequence']
            ids["inchikey"] = row['standard_inchi_key']
            ids["chembl_ligand"] = row['molecule_chembl_id']

            # Filter types
            ligand = get_or_create_ligand(row['pref_name'], ids, ligand_types[row['molecule_type']], False, False)

            # Add LigandIDs
            if row['other_ids'] != None:
                extra_ids = row['other_ids'].split(";")
                existing = list(set.intersection(set(extra_ids), set(existing_ids)))
                if len(existing) > 0:
                    continue # skip rest of creation
                else:
                    for id in extra_ids:
                        weblinks.append(LigandID(ligand = ligand, index=id, web_resource=wr_chembl))
        # Bulk insert all new ligandIDs
        LigandID.objects.bulk_create(weblinks)


    @staticmethod
    def build_chembl_bioactivities():
        AssayExperiment.objects.all().delete()
        print("\n===============\n#1 Reading ChEMBL bioacitivity data")
        bioactivity_input_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_bioactivity_data.csv"])
        bioactivity_data = pd.read_csv(bioactivity_input_file, dtype=str)
        bio_entries = len(bioactivity_data)
        print("Found", bio_entries, "bioactivities", datetime.datetime.now())

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        #ids = list(bioactivity_data["parent_molecule_chembl_id"].unique())  # not filtering is way faster
        ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Building ChEMBL proteins cache", datetime.datetime.now())
        # NOTE => might need to switch to Accession as the Entry name changes more frequently
        # If so, keep isoform notations in mind
        names = list(bioactivity_data["Entry name"].unique())
        proteins = list(Protein.objects.filter(entry_name__in=names).values_list("pk", "entry_name"))
        prot_dict = {prot_entry[1]: prot_entry[0] for prot_entry in proteins}

        print("\n#4 Building ChEMBL bioactivity entries", datetime.datetime.now())
        bioacts = []
        for index, row in bioactivity_data.iterrows():
            bioacts.append(AssayExperiment())
            bioacts[-1].ligand_id = lig_dict[row["parent_molecule_chembl_id"]]
            bioacts[-1].protein_id = prot_dict[row["Entry name"]]
            bioacts[-1].assay_type = row["assay_type"]
            bioacts[-1].assay_description = row["assay_description"]
            bioacts[-1].pchembl_value = row["pchembl_value"]
            bioacts[-1].standard_value = row["standard_value"]
            bioacts[-1].standard_relation = row["standfard_relation"]
            bioacts[-1].standard_type = row["standard_type"]
            bioacts[-1].standard_units = row["standard_units"]
            bioacts[-1].document_chembl_id = row["document_chembl_id"]

            # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of", bio_entries, "bioactivities")
                bioacts = []


    @staticmethod
    def build_pubchem_vendor_links():
        LigandVendors.objects.all().delete()
        print("\n===============\n#1 Reading and creating Vendors")
        vendor_url = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_list.csv"])
        vendor_data = pd.read_csv(vendor_url, dtype=str)
        vendors = []
        for index, row in vendor_data.iterrows():
            vendors.append(LigandVendors(slug=slugify(row["SourceName"]), name = row["SourceName"], url = row["SourceURL"]))
        LigandVendors.objects.bulk_create(vendors)
        vendor_dict = {vendor.name : vendor.pk for vendor in vendors}

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Creating all vendor links", datetime.datetime.now())
        vendor_links_url = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_links.csv"])
        vendor_links_data = pd.read_csv(vendor_links_url, dtype=str)
        links = []
        for index, row in vendor_links_data.iterrows():
            if len(row["SourceRecordURL"]) < 300:
                links.append(LigandVendorLink(vendor_id=vendor_dict[row["SourceName"]], ligand_id = lig_dict[row["chembl_id"]], url = row["SourceRecordURL"]))

        LigandVendorLink.objects.bulk_create(links)
