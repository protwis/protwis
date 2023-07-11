from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from structure.models import Structure
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd


import logging, json, os
from collections import OrderedDict



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)

    pdbs = '''    6CM4    Structure of the D2 Dopamine Receptor Bound to the Atypical Antipsychotic Drug Risperidone  3/14/2018
    6FFH    Crystal Structure of mGluR5 in complex with Fenobam at 2.65 A   3/7/2018
    6FFI    Crystal Structure of mGluR5 in complex with MMPEP at 2.2 A  3/7/2018
    5WF5    Agonist bound human A2a adenosine receptor with D52N mutation at 2.60 A resolution  2/21/2018
    5WF6    Agonist bound human A2a adenosine receptor with S91A mutation at 2.90 A resolution  2/21/2018
    6BQG    Crystal structure of 5-HT2C in complex with ergotamine  2/14/2018
    6BQH    Crystal structure of 5-HT2C in complex with ritanserin  2/14/2018
    5V54    Crystal structure of 5-HT1B receptor in complex with methiothepin   2/7/2018
    5OLG    Structure of the A2A-StaR2-bRIL562-ZM241385 complex at 1.86A obtained from in meso soaking experiments. 1/17/2018
    5OLH    Structure of the A2A-StaR2-bRIL562-Vipadenant complex at 2.6A obtained from in meso soaking experiments.    1/17/2018
    5OLO    Structure of the A2A-StaR2-bRIL562-Tozadenant complex at 3.1A obtained from in meso soaking experiments.    1/17/2018
    5OLV    Structure of the A2A-StaR2-bRIL562-LUAA47070 complex at 2.0A obtained from in meso soaking experiments. 1/17/2018
    5OLZ    Structure of the A2A-StaR2-bRIL562-Compound 4e complex at 1.9A obtained from bespoke co-crystallisation experiments.    1/17/2018
    5OM1    Structure of the A2A-StaR2-bRIL562-Compound 4e complex at 2.1A obtained from in meso soaking experiments (1 hour soak). 1/17/2018
    5OM4    Structure of the A2A-StaR2-bRIL562-Compound 4e complex at 1.86A obtained from in meso soaking experiments (24 hour soak).   1/17/2018
    6B73    Crystal Structure of a nanobody-stabilized active state of the kappa-opioid receptor    1/17/2018
    5O9H    Crystal structure of thermostabilised human C5a anaphylatoxin chemotactic receptor 1 (C5aR) in complex with NDT9513727  1/10/2018
    6AQF    Crystal structure of A2AAR-BRIL in complex with the antagonist ZM241385 produced from Pichia pastoris   1/10/2018
    5VRA    2.35-Angstrom In situ Mylar structure of human A2A adenosine receptor at 100 K  12/13/2017
    5WKT    3.2-Angstrom In situ Mylar structure of bovine opsin at 100 K   12/13/2017
    5WIU    Structure of the human D4 Dopamine receptor in complex with Nemonapride 10/18/2017
    5WIV    Structure of the sodium-bound human D4 Dopamine receptor in complex with Nemonapride    10/18/2017
    5NLX    A2A Adenosine receptor room-temperature structure determined by serial millisecond crystallography  9/27/2017
    5X7D    Structure of beta2 adrenoceptor bound to carazolol and an intracellular allosteric antagonist   8/16/2017
    5W0P    Crystal structure of rhodopsin bound to visual arrestin determined by X-ray free electron laser 8/9/2017
    5TUD    Structural Insights into the Extracellular Recognition of the Human Serotonin 2B Receptor by an Antibody    7/26/2017
    5NX2    Crystal structure of thermostabilised full-length GLP-1R in complex with a truncated peptide agonist at 3.7 A resolution    6/14/2017
    5VBL    Structure of apelin receptor in complex with agonist peptide    5/31/2017
    5UVI    Serial Millisecond Crystallography of Membrane and Soluble Protein Micro-crystals using Synchrotron Radiation   5/24/2017
    5V56    2.9A XFEL structure of the multi-domain human smoothened receptor (with E194M mutation) in complex with TC114   5/24/2017
    5V57    3.0A SYN structure of the multi-domain human smoothened receptor in complex with TC114  5/24/2017
    5VAI    Cryo-EM structure of the activated Glucagon-like peptide-1 receptor in complex with G protein   5/24/2017
    5VEW    Structure of the human GLP-1 receptor complex with PF-06372222  5/24/2017
    5VEX    Structure of the human GLP-1 receptor complex with NNC0640  5/17/2017
    5NDD    Crystal structure of a thermostabilised human protease-activated receptor-2 (PAR2) in complex with AZ8838 at 2.8 angstrom resolution    5/3/2017
    5NDZ    Crystal structure of a thermostabilised human protease-activated receptor-2 (PAR2) in complex with AZ3451 at 3.6 angstrom resolution    5/3/2017
    5NJ6    Crystal structure of a thermostabilised human protease-activated receptor-2 (PAR2) in ternary complex with Fab3949 and AZ7188 at 4.0 angstrom resolution    5/3/2017
    5UZ7    Volta phase plate cryo-electron microscopy structure of a calcitonin receptor-heterotrimeric Gs protein complex 5/3/2017
    5UNF    XFEL structure of human angiotensin II type 2 receptor (Monoclinic form) in complex with compound 1 (N-benzyl-N-(2-ethyl-4-oxo-3-{[2'-(2H-tetrazol-5-yl)[1,1'-biphenyl]-4-yl])  4/5/2017
    5UNG    XFEL structure of human angiotensin II type 2 receptor (Orthorhombic form) in complex with compound 1 (N-benzyl-N-(2-ethyl-4-oxo-3-{[2'-(2H-tetrazol-5-yl)[1,1'-biphenyl]-4-yl] methyl}-3,4-dihydroquinazolin-6-yl)thiophene-2-carboxamide) 4/5/2017
    5UNH    Synchrotron structure of human angiotensin II type 2 receptor in complex with compound 2 (N-[(furan-2-yl)methyl]-N-(4-oxo-2-propyl-3-{[2'-(2H-tetrazol-5-yl)[1,1'- biphenyl]-4-yl]methyl}-3,4-dihydroquinazolin-6-yl)benzamide) 4/5/2017
    5UEN    Crystal structure of the human adenosine A1 receptor A1AR-bRIL in complex with the covalent antagonist DU172 at 3.2A resolution 3/1/2017
    5UIG    Crystal structure of adenosine A2A receptor bound to a novel triazole-carboximidamide antagonist    2/8/2017
    5TVN    Crystal structure of the LSD-bound 5-HT2B receptor  2/1/2017
    5T04    STRUCTURE OF CONSTITUTIVELY ACTIVE NEUROTENSIN RECEPTOR 12/21/2016
    5T1A    Structure of CC Chemokine Receptor 2 with Orthosteric and Allosteric Antagonists    12/14/2016
    5TDH    The crystal structure of the dominant negative mutant G protein alpha(i)-1-beta-1-gamma-2 G203A/A326S   11/9/2016
    5KVM    Extracellular region of mouse GPR56/ADGRG1 in complex with FN3 monobody 9/28/2016
    5K2A    2.5 angstrom A2a adenosine receptor structure with sulfur SAD phasing using XFEL data   9/21/2016
    5K2B    2.5 angstrom A2a adenosine receptor structure with MR phasing using XFEL data   9/21/2016
    5K2C    1.9 angstrom A2a adenosine receptor structure with sulfur SAD phasing and phase extension using XFEL data   9/21/2016
    5K2D    1.9A angstrom A2a adenosine receptor structure with MR phasing using XFEL data  9/21/2016
    5D6L    beta2AR-T4L - CIM   8/17/2016
    5KZV    Crystal structure of the xenopus Smoothened cysteine-rich domain (CRD) in complex with 20(S)-hydroxycholesterol 8/17/2016
    5KZY    Crystal structure of the xenopus Smoothened cysteine-rich domain (CRD) in complex with cyclopamine  8/17/2016
    5KZZ    Crystal structure of the xenopus Smoothened cysteine-rich domain (CRD) in its apo-form  8/17/2016
    5G53    Structure of the adenosine A2A receptor bound to an engineered G protein    8/3/2016
    5JQH    Structure of beta2 adrenoceptor bound to carazolol and inactive-state stabilizing nanobody, Nb60    7/13/2016
    4Z9G    Crystal structure of human corticotropin-releasing factor receptor 1 (CRF1R) in complex with the antagonist CP-376395 in a hexagonal setting with translational non-crystallographic symmetry   6/29/2016
    5JS7    Structural model of a apo G-protein alpha subunit determined with NMR residual dipolar couplings and SAXS   6/29/2016
    5JS8    Structural Model of a Protein alpha subunit in complex with GDP obtained with SAXS and NMR residual couplings   6/29/2016
    5FTT    Octameric complex of Latrophilin 3 (Lec, Olf) , Unc5D (Ig, Ig2, TSP1) and FLRT2 (LRR)   5/4/2016
    5FTU    Tetrameric complex of Latrophilin 3, Unc5D and FLRT2    5/4/2016
    2N55    Structure of constitutively monomeric CXCL12 in complex with the CXCR4 N-terminus   4/27/2016
    5EE7    Crystal structure of the human glucagon receptor (GCGR) in complex with the antagonist MK-0893  4/20/2016
    5DGY    Crystal structure of rhodopsin bound to visual arrestin 3/23/2016
    5DSG    Structure of the M4 muscarinic acetylcholine receptor (M4-mT4L) bound to tiotropium 3/16/2016
    5CXV    Structure of the human M1 muscarinic acetylcholine receptor bound to antagonist Tiotropium  3/9/2016
    4ZRG    Visual arrestin mutant - R175E  11/11/2015
    4X1H    Opsin/G(alpha) peptide complex stabilized by nonyl-glucoside    11/4/2015
    5DHG    The crystal structure of nociceptin/orphanin FQ peptide receptor (NOP) in complex with C-35 (PSI Community Target)  10/21/2015
    5DHH    The crystal structure of nociceptin/orphanin FQ peptide receptor (NOP) in complex with SB-612111 (PSI Community Target) 10/21/2015
    4ZUD    Crystal Structure of Human Angiotensin Receptor in Complex with Inverse Agonist Olmesartan at 2.8A resolution.  10/7/2015
    2N2F    Solution NMR structure of Dynorphin 1-13 bound to Kappa Opioid Receptor 9/9/2015
    4RMK    Crystal structure of the Olfactomedin domain of latrophilin 3 in P65 crystal form   8/19/2015
    4RML    Crystal structure of the Olfactomedin domain of latrophilin 3 in C2221 crystal form 8/19/2015
    5CGC    Structure of the human class C GPCR metabotropic glutamate receptor 5 transmembrane domain in complex with the negative allosteric modulator 3-chloro-4-fluoro-5-[6-(1H-pyrazol-1-yl)pyrimidin-4-yl]benzonitrile    8/12/2015
    5CGD    Structure of the human class C GPCR metabotropic glutamate receptor 5 transmembrane domain in complex with the negative allosteric modulator 3-chloro-5-[6-(5-fluoropyridin-2-yl)pyrimidin-4-yl]benzonitrile - (HTL14242)   8/12/2015
    4XEE    Structure of active-like neurotensin receptor   7/29/2015
    4XES    Structure of active-like neurotensin receptor   7/29/2015
    4ZWJ    Crystal structure of rhodopsin bound to arrestin by femtosecond X-ray laser 7/29/2015
    4TNB    Crystal Structure of G Protein-Coupled Receptor Kinase 5 in Complex with Sangivamycin   6/10/2015
    4TND    Crystal Structure of G Protein-Coupled Receptor Kinase 5 in Complex with AMP-PNP    6/10/2015
    4Z34    Crystal Structure of Human Lysophosphatidic Acid Receptor 1 in complex with ONO9780307  6/3/2015
    4Z35    Crystal Structure of Human Lysophosphatidic Acid Receptor 1 in complex with ONO-9910539 6/3/2015
    4Z36    Crystal Structure of Human Lysophosphatidic Acid Receptor 1 in complex with ONO-3080573 6/3/2015
    4YAY    XFEL structure of human Angiotensin Receptor    4/22/2015
    4UG2    Thermostabilised HUMAN A2a Receptor with CGS21680 bound 4/8/2015
    4UHR    Thermostabilised HUMAN A2a Receptor with CGS21680 bound 4/8/2015
    4XNV    The human P2Y1 receptor in complex with BPTU    4/1/2015
    4XNW    The human P2Y1 receptor in complex with MRS2500 4/1/2015
    4XT1    Structure of a nanobody-bound viral GPCR bound to human chemokine CX3CL1    3/4/2015
    4XT3    Structure of a viral GPCR bound to human chemokine CX3CL1   3/4/2015
    4RWS    Crystal structure of CXCR4 and viral chemokine antagonist vMIP-II complex (PSI Community Target)    2/11/2015
    4RWA    Synchrotron structure of the human delta opioid receptor in complex with a bifunctional peptide (PSI community target)  1/14/2015
    4RWD    XFEL structure of the human delta opioid receptor in complex with a bifunctional peptide    1/14/2015
    4UV4    Crystal structure of anti-FPR Fpro0165 Fab fragment 12/24/2014
    4U14    Structure of the M3 muscarinic acetylcholine receptor bound to the antagonist tiotropium crystallized with disulfide-stabilized T4 lysozyme (dsT4L) 11/26/2014
    4U15    M3-mT4L receptor bound to tiotropium    11/26/2014
    4U16    M3-mT4L receptor bound to NMS   11/26/2014
    4PNI    Bovine G protein-coupled receptor kinase 1 in complex with GSK2163632A  10/8/2014
    4R7V    Crystal structure of N-lobe of human ARRDC3(1-165)  10/1/2014
    4R7X    Crystal structure of N-lobe of human ARRDC3(1-180)  10/1/2014
    4PXF    Crystal structure of the active G-protein-coupled receptor opsin in complex with the finger-loop peptide derived from the full-length arrestin-1    9/17/2014
    4QIM    Structure of the human smoothened receptor in complex with ANTA XV  7/23/2014
    4QIN    Structure of the human smoothened receptor in complex with SAG1.5   7/23/2014
    4OO9    Structure of the human class C GPCR metabotropic glutamate receptor 5 transmembrane domain in complex with the negative allosteric modulator mavoglurant    7/2/2014
    4P39    Crystal structure of the human C5aR antagonist C5a-A8   6/11/2014
    4P3A    Crystal structure of the mouse C5a anaphylatoxin    6/11/2014
    4P3B    Crystal structure of the mouse C5a-desArg anaphylatoxin 6/11/2014
    4PXZ    Crystal structure of P2Y12 receptor in complex with 2MeSADP 4/30/2014
    4PY0    Crystal structure of P2Y12 receptor in complex with 2MeSATP 4/30/2014
    4MQW    Structure of follicle-stimulating hormone in complex with the entire ectodomain of its receptor (P31)   4/9/2014
    4NTJ    Structure of the human P2Y12 receptor in complex with an antithrombotic drug    3/26/2014
    4OR2    Human class C G protein-coupled metabotropic glutamate receptor 1 in complex with a negative allosteric modulator   3/19/2014
    4O9R    Human Smoothened Receptor structure in complex with cyclopamine 3/5/2014
    4NUU    Heterotrimer structure of Region II from Plasmodium vivax Duffy Binding Protein (PvDBP) bound to the ectodomain of the Duffy Antigen Receptor for Chemokines (DARC) 2/5/2014
    4NUV    Heterotetramer structure of Region II from Plasmodium vivax Duffy Binding Protein (PvDBP) bound to the ectodomain of the Duffy Antigen Receptor for Chemokines (DARC)   2/5/2014
    4L9I    Bovine G Protein Coupled Receptor Kinase 1 in Complex with Paroxetine   1/22/2014
    4N4W    Structure of the human smoothened receptor in complex with SANT-1.  1/22/2014
    4N6H    1.8 A Structure of the human delta opioid 7TM receptor (PSI Community Target)   12/25/2013
    4NC3    Crystal structure of the 5-HT2B receptor solved using serial femtosecond crystallography in lipidic cubic phase.    12/18/2013
    2M9O    Solution structure of kalata B7 11/20/2013
    4J4Q    Crystal structure of active conformation of GPCR opsin stabilized by octylglucoside 10/30/2013
    4MBS    Crystal Structure of the CCR5 Chemokine Receptor    9/11/2013
    4L6R    Structure of the class B human glucagon G protein coupled receptor  7/24/2013
    4K5Y    Crystal structure of human corticotropin-releasing factor receptor 1 (CRF1R) in complex with the antagonist CP-376395   7/17/2013
    4BEY    Night blindness causing G90D rhodopsin in complex with GaCT2 peptide    5/8/2013
    4J2Q    Crystal structure of C-terminally truncated arrestin reveals mechanism of arrestin activation   5/1/2013
    4BEZ    Night blindness causing G90D rhodopsin in the active conformation   4/24/2013
    4JKV    Structure of the human smoothened 7TM receptor in complex with an antitumor agent   4/24/2013
    2LMK    Solution Structure of Mouse Pheromone ESP1  4/17/2013
    3T33    Crystal Structure of Arabidopsis GCR2   4/17/2013
    4JQI    Structure of active beta-arrestin1 bound to a G protein-coupled receptor phosphopeptide 4/17/2013
    3ZPQ    Thermostabilised turkey beta1 adrenergic receptor with 4-(piperazin-1- yl)-1H-indole bound (compound 19)    4/3/2013
    3ZPR    Thermostabilised turkey beta1 adrenergic receptor with 4-methyl-2-(piperazin-1-yl) quinoline bound  4/3/2013
    4IAQ    Crystal structure of the chimeric protein of 5-HT1B-BRIL in complex with dihydroergotamine (PSI Community Target)   3/13/2013
    4IAR    Crystal structure of the chimeric protein of 5-HT1B-BRIL in complex with ergotamine (PSI Community Target)  3/13/2013
    4IB4    Crystal structure of the chimeric protein of 5-HT2B-BRIL in complex with ergotamine 3/13/2013
    4I6O    Crystal structure of chemically synthesized human anaphylatoxin C3a 2/27/2013
    2LQ4    Structural Characterization of an LPA1 Second Extracellular Loop Mimetic with a Self-Assembling Coiled-Coil Folding Constraint  1/9/2013
    2LNL    Structure of human CXCR1 in phospholipid bilayers   10/17/2012
    4ERS    A Molecular Basis for Negative Regulation of the Glucagon Receptor  8/29/2012
    4AY9    Structure of follicle-stimulating hormone in complex with the entire ectodomain of its receptor 8/8/2012
    4EIY    Crystal structure of the chimeric protein of A2aAR-BRIL in complex with ZM241385 at 1.8A resolution 7/25/2012
    3T8O    Rhodopsin kinase (GRK1) L166K mutant at 2.5A resolution 6/6/2012
    4EA3    Structure of the N/OFQ Opioid Receptor in Complex with a Peptide Mimetic    4/25/2012
    2RRS    NMR Structure of LC4 transmembrane segment of CCR5  4/11/2012
    3UZA    Thermostabilised Adenosine A2A receptor in complex with 6-(2,6-Dimethylpyridin-4-yl)-5-phenyl-1,2,4-triazin-3-amine 3/21/2012
    3UZC    Thermostabilised Adenosine A2A receptor in complex with 4-(3-amino-5-phenyl-1,2,4-triazin-6-yl)-2-chlorophenol  3/21/2012
    4DJH    Structure of the human kappa opioid receptor in complex with JDTic  3/21/2012
    4DLO    Crystal structure of the GAIN and HormR domains of brain angiogenesis inhibitor 3 (BAI3)    2/22/2012
    4DLQ    Crystal structure of the GAIN and HormR domains of CIRL 1/Latrophilin 1 (CL1)   2/22/2012
    3V2W    Crystal Structure of a Lipid G protein-Coupled Receptor at 3.35A    2/15/2012
    3V2Y    Crystal Structure of a Lipid G protein-Coupled Receptor at 2.80A    2/15/2012
    3UGU    Crystal Structure of p44 (Splice Variant of Visual Arrestin)    2/8/2012
    3UGX    Crystal Structure of Visual Arrestin    2/8/2012
    3UMS    Crystal structure of the G202A mutant of human G-alpha-i1   2/8/2012
    3UON    Structure of the human M2 muscarinic acetylcholine receptor bound to an antagonist  2/1/2012
    3AQE    Crystal structure of the extracellular domain of human RAMP2    11/9/2011
    3AQF    Crystal structure of the human CRLR/RAMP2 extracellular complex 11/2/2011
    2L60    A novel design concept: New Y-receptor agonists with increased membrane recruitment, Y2 affinity and selectivity    10/5/2011
    3PWH    Thermostabilised Adenosine A2A Receptor 9/7/2011
    3REY    Thermostabilised adenosine A2A receptor in complex with XAC 9/7/2011
    3RFM    Thermostabilised adenosine A2A receptor in complex with caffeine    9/7/2011
    3SN6    Crystal structure of the beta2 adrenergic receptor-Gs protein complex   7/20/2011
    3RZE    Structure of the human histamine H1 receptor in complex with doxepin    6/15/2011
    2YCY    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND ANTAGONIST CYANOPINDOLOL  6/8/2011
    2YCW    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND ANTAGONIST CARAZOLOL  6/1/2011
    2YCX    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND ANTAGONIST CYANOPINDOLOL  6/1/2011
    2YCZ    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND ANTAGONIST IODOCYANOPINDOLOL  6/1/2011
    2YDO    Thermostabilised HUMAN A2a Receptor with adenosine bound    5/18/2011
    2YDV    Thermostabilised HUMAN A2a Receptor with NECA bound 5/18/2011
    2Y01    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND PARTIAL AGONIST DOBUTAMINE (CRYSTAL DOB102)   3/30/2011
    2XWT    CRYSTAL STRUCTURE OF THE TSH RECEPTOR IN COMPLEX WITH A BLOCKING TYPE TSHR AUTOANTIBODY 3/9/2011
    3QAK    Agonist bound structure of the human adenosine A2a receptor 3/9/2011
    2L63    NMR solution structure of GLP-2 in 2,2,2 trifluroethanol    2/23/2011
    2L64    NMR Solution structure of GLP-2 in DHPC micelles    2/23/2011
    3NKM    Crystal structure of mouse autotaxin    1/19/2011
    3NKN    Crystal structure of mouse autotaxin in complex with 14:0-LPA   1/19/2011
    3NKO    Crystal structure of mouse autotaxin in complex with 16:0-LPA   1/19/2011
    3NKP    Crystal structure of mouse autotaxin in complex with 18:1-LPA   1/19/2011
    3NKQ    Crystal structure of mouse autotaxin in complex with 18:3-LPA   1/19/2011
    3NKR    Crystal structure of mouse autotaxin in complex with 22:6-LPA   1/19/2011
    3P0G    Structure of a nanobody-stabilized active state of the beta2 adrenoceptor   1/19/2011
    2Y00    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND PARTIAL AGONIST DOBUTAMINE (CRYSTAL DOB92)    1/12/2011
    2Y02    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND AGONIST CARMOTEROL    1/12/2011
    2Y03    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND AGONIST ISOPRENALINE  1/12/2011
    2Y04    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND PARTIAL AGONIST SALBUTAMOL    1/12/2011
    3PDS    Irreversible Agonist-Beta2 Adrenoceptor Complex 1/12/2011
    2XVT    Structure of the extracellular domain of human RAMP2    12/29/2010
    3PBL    Structure of the human dopamine D3 receptor in complex with eticlopride 11/3/2010
    3ODU    The 2.5 A structure of the CXCR4 chemokine receptor in complex with small molecule antagonist IT1t  10/27/2010
    3OE0    Crystal structure of the CXCR4 chemokine receptor in complex with a cyclic peptide antagonist CVX15 10/27/2010
    3OE6    Crystal structure of the CXCR4 chemokine receptor in complex with a small molecule antagonist IT1t in I222 spacegroup   10/27/2010
    3OE8    Crystal structure of the CXCR4 chemokine receptor in complex with a small molecule antagonist IT1t in P1 spacegroup 10/27/2010
    3OE9    Crystal structure of the chemokine CXCR4 receptor in complex with a small molecule antagonist IT1t in P1 spacegroup 10/27/2010
    3N93    Crystal structure of human CRFR2 alpha extracellular domain in complex with Urocortin 3 10/20/2010
    3N95    Crystal structure of human CRFR2 alpha extracellular domain in complex with Urocortin 2 10/20/2010
    3N96    Crystal structure of human CRFR2 alpha extracellular domain in complex with Urocortin 1 10/20/2010
    3N7P    Crystal structure of the ectodomain complex of the CGRP receptor, a Class-B GPCR, reveals the site of drug antagonism   9/15/2010
    3N7R    Crystal structure of the ectodomain complex of the CGRP receptor, a Class-B GPCR, reveals the site of drug antagonism   9/15/2010
    3N7S    Crystal structure of the ectodomain complex of the CGRP receptor, a Class-B GPCR, reveals the site of drug antagonism   9/15/2010
    3NY8    Crystal structure of the human beta2 adrenergic receptor in complex with the inverse agonist ICI 118,551    8/11/2010
    3NY9    Crystal structure of the human beta2 adrenergic receptor in complex with a novel inverse agonist    8/11/2010
    3NYA    Crystal structure of the human beta2 adrenergic receptor in complex with the neutral antagonist alprenolol  8/11/2010
    2X57    Crystal structure of the extracellular domain of human Vasoactive intestinal polypeptide receptor 2 3/9/2010
    2KOE    Human cannabinoid receptor 1 - helix 7/8 peptide    10/6/2009
    3H3G    Crystal structure of the extracellular domain of the human parathyroid hormone receptor (PTH1R) in complex with parathyroid hormone-related protein (PTHrP) 8/11/2009
    3G04    Crystal structure of the TSH receptor in complex with a thyroid-stimulating autoantibody    8/4/2009
    3I8S    Structure of the cytosolic domain of E. coli FeoB, nucleotide-free form 7/28/2009
    3I8X    Structure of the cytosolic domain of E. coli FeoB, GDP-bound form   7/28/2009
    3I92    Structure of the cytosolic domain of E. coli FeoB, GppCH2p-bound form   7/28/2009
    2KI9    Human cannabinoid receptor-2 helix 6    5/12/2009
    2K9P    Structure of TM1_TM2 in LPPG micelles   5/5/2009
    2ZJY    Structure of the K349P mutant of Gi alpha 1 subunit bound to ALF4 and GDP   3/24/2009
    2ZJZ    Structure of the K349P mutant of Gi alpha 1 subunit bound to GDP    3/24/2009
    2K3U    Structure of the tyrosine-sulfated C5a receptor N-terminus in complex with the immune evasion protein CHIPS.    3/10/2009
    2JX4    NMR structure of the intracellular loop (i3) of the vasopressin V2 receptor (GPCR)  11/18/2008
    3EML    The 2.6 A Crystal Structure of a Human A2A Adenosine Receptor bound to ZM241385.    10/14/2008
    2VT4    TURKEY BETA1 ADRENERGIC RECEPTOR WITH STABILISING MUTATIONS AND BOUND CYANOPINDOLOL 6/24/2008
    3D4S    Cholesterol bound form of human beta2 adrenergic receptor.  6/17/2008
    2JX0    The paxillin-binding domain (PBD) of G Protein Coupled Receptor (GPCR)-kinase (GRK) interacting protein 1 (GIT1)    4/29/2008
    2JX9    Solution structure of the Gal_lectin domain of mouse Latrophilin-1 GPCR 4/29/2008
    2JXA    Mouse Latrophilin-1 GPCR Gal_lectin domain in complex with Rhamnose 4/29/2008
    2YX8    Crystal structure of the extracellular domain of human RAMP1    4/29/2008
    2RGN    Crystal Structure of p63RhoGEF complex with Galpha-q and RhoA   1/15/2008
    2NYZ    Viral Chemokine Binding Protein M3 From Murine Gammaherpesvirus68 In Complex With The C- Chemokine XCL1 12/25/2007
    2NZ1    Viral Chemokine Binding Protein M3 From Murine Gammaherpesvirus68 In Complex With The CC-Chemokine CCL2/MCP-1   12/25/2007
    2RH1    High resolution crystal structure of human B2-adrenergic G protein-coupled receptor.    10/30/2007
    2PZX    Structure of the methuselah ectodomain with peptide inhibitor   8/28/2007
    2QKH    Crystal structure of the extracellular domain of human GIP receptor in complex with the hormone GIP 8/14/2007
    2DCO    S1P4 First Extracellular Loop Peptidomimetic    1/23/2007
    2O8Z    Bound Structure of CRF1 Extracellular Domain Antagonist 12/26/2006
    2I37    Crystal structure of a photoactivated rhodopsin 10/17/2006
    1YM7    G Protein-Coupled Receptor Kinase 2 (GRK2)  7/5/2005
    1YTV    Maltose-binding protein fusion to a C-terminal fragment of the V1a vasopressin receptor 4/12/2005
    1Y7J    NMR structure family of Human Agouti Signalling Protein (80-132: Q115Y, S124Y)  2/15/2005
    1Y7K    NMR structure family of Human Agouti Signalling Protein (80-132: Q115Y, S124Y)  2/15/2005
    1WSO    The solution structures of human Orexin-A   11/30/2004
    1NZS    NMR structures of phosphorylated carboxy terminus of bovine rhodopsin in arrestin-bound state   3/2/2004
    1MF6    Transducin gamma subunit, C-terminal domain 60-71, rhodopsin-bound state: Ensemble of 15 models determined by TrNOE spectroscopy    4/15/2003
    1HZX    CRYSTAL STRUCTURE OF BOVINE RHODOPSIN   7/4/2001
    1FJR    CRYSTAL STRUCTURE OF THE ECTODOMAIN OF METHUSELAH   4/4/2001
    1D6G    MOLECULAR COMPLEX OF CHOLECYSTOKININ-8 AND N-TERMINUS OF THE CHOLECYSTOKININ A RECEPTOR BY NMR SPECTROSCOPY 11/17/1999'''


    def handle(self, *args, **options):

        # for line in self.pdbs.splitlines():
        #     t = line.split("    ")
        #     # print(t[1])
        #     exists = Structure.objects.filter(pdb_code__index = t[1]).exists()
        #     if not exists:
        #         print(t[1],"not there",t)
        proteins = set()
        proteins_nonortho = set()
        p2pdb = {}
        pdb2p = {}
        new_pdbs = ['5O9H', '5OLG', '5OLH', '5OLO', '5OLV', '5OLZ', '5OM1', '5OM4', '5V54', '5WF5', '5WF6', '5YQZ', '6AQF', '6B3J', '6B73', '6BQG', '6BQH', '6CM4', '6FFH', '6FFI']
        for s in Structure.objects.all().exclude(structure_type__slug__startswith='af-'):
           # print(s)
            #print(s.protein_conformation.protein.parent.entry_name)
            proteins.add(s.protein_conformation.protein.parent.entry_name)
            proteins_nonortho.add(s.protein_conformation.protein.parent.family.name)
            if s.protein_conformation.protein.parent.entry_name not in p2pdb:
                p2pdb[s.protein_conformation.protein.parent.entry_name] = []
            p2pdb[s.protein_conformation.protein.parent.entry_name].append(str(s))
            pdb2p[str(s)] = s.protein_conformation.protein.parent.entry_name
        print(len(proteins),"total entry names")
        print(len(proteins_nonortho),"total entry names (proteins_nonortho)")

        for p in new_pdbs:
            print('New PDB',p,'entry_name',pdb2p[p],'other pdbs for this protein',p2pdb[pdb2p[p]])
