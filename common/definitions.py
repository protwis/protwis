'''Common definitions to be used throughout the Django module.'''

from collections import OrderedDict
from decimal import Decimal

AMINO_ACIDS = OrderedDict([
    ('A', 'Ala'),
    ('C', 'Cys'),
    ('D', 'Asp'),
    ('E', 'Glu'),
    ('F', 'Phe'),
    ('G', 'Gly'),
    ('H', 'His'),
    ('I', 'Ile'),
    ('K', 'Lys'),
    ('L', 'Leu'),
    ('M', 'Met'),
    ('N', 'Asn'),
    ('P', 'Pro'),
    ('Q', 'Gln'),
    ('R', 'Arg'),
    ('S', 'Ser'),
    ('T', 'Thr'),
    ('V', 'Val'),
    ('W', 'Trp'),
    ('Y', 'Tyr'),
    ('B', 'Asx'),
    ('Z', 'Glx'),
    ('J', 'Xle')
])

# Amino Acid Propensities
AA_PROPENSITY = {
    "A":Decimal('0'),
    "R":Decimal('0.21'),
    "N":Decimal('0.65'),
    "D":Decimal('0.69'),
    "C":Decimal('0.68'),
    "E":Decimal('0.40'),
    "Q":Decimal('0.39'),
    "G":Decimal('1.00'),
    "H":Decimal('0.61'),
    "I":Decimal('0.41'),
    "L":Decimal('0.21'),
    "K":Decimal('0.26'),
    "M":Decimal('0.24'),
    "F":Decimal('0.54'),
    "P":Decimal('3.16'),
    "S":Decimal('0.50'),
    "T":Decimal('0.66'),
    "W":Decimal('0.49'),
    "Y":Decimal('0.53'),
    "V":Decimal('0.61')
}

# The Octanol-Interface scale from the Wimley-White scale uis used here.
HYDROPHOBICITY = {
    "I":Decimal('-1.12'),
    "L":Decimal('-1.25'),
    "F":Decimal('-1.71'),
    "V":Decimal('-0.46'),
    "M":Decimal('-0.67'),
    "P":Decimal('0.14'),
    "W":Decimal('-2.09'),
    "H":Decimal('0.11'),
    "T":Decimal('0.25'),
    "E":Decimal('0.11'),
    "Q":Decimal('0.77'),
    "C":Decimal('-0.02'),
    "Y":Decimal('-0.71'),
    "S":Decimal('0.46'),
    "A":Decimal('0.50'),
    "N":Decimal('0.85'),
    "D":Decimal('0.43'),
    "R":Decimal('1.81'),
    "G":Decimal('1.15'),
    "K":Decimal('2.80'),
    "R+":Decimal('1.81'),
    "H+":Decimal('2.33'),
    "E-":Decimal('3.63'),
    "K+":Decimal('2.80'),
    "D-":Decimal('3.64')
    }


AMINO_ACID_GROUPS = OrderedDict([
        ('hp',     ('A', 'C', 'F', 'I', 'L', 'M', 'P', 'V', 'W', 'Y')),
        ('alhp',   ('A', 'C', 'I', 'L', 'M', 'P', 'V')),
        ('arhp',   ('F', 'W', 'Y')),
        ('ar',     ('F', 'H', 'W', 'Y')),
        ('pol',    ('D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T')),
        ('hb',     ('D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y')),
        ('hbd',    ('H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y')),
        ('hba',    ('D', 'E', 'H', 'N', 'Q', 'S', 'T', 'Y')),
        ('charge', ('D', 'E', 'H', 'K', 'R')),
        ('neg',    ('D', 'E')),
        ('pos',    ('H', 'K', 'R')),
        ('lar',    ('E', 'F', 'H', 'K', 'Q', 'R', 'W', 'Y')),
        ('sma',    ('A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V')),
        ('any',    ()),
        ('custom', ()),
    ])

AMINO_ACID_GROUP_NAMES = OrderedDict([
        ('hp',     'Hydrophobic - All'),
        ('alhp',   'Hydrophobic - Aliphatic'),
        ('arhp',   'Hydrophobic - Aromatic'),
        ('ar',     'Aromatic'),
        ('pol',    'Polar'),
        ('hb',     'H-Bonding'),
        ('hbd',    'H-Bond Donor'),
        ('hba',    'H-Bond Acceptor'),
        ('charge', 'Charge'),
        ('neg',    'Negative charge'),
        ('pos',    'Positive charge'),
        ('lar',    'Large'),
        ('sma',    'Small'),
        ('any',    'Any feature'),
        ('custom', 'Custom'),
    ])

DESIGN_SUBSTITUTION_MATRIX = OrderedDict([
        ('hydrophobic', {
            #'A':[[['L','V','I']],['Increase size to "block" binding']],
            'C':[['A'],['Remove vdW interactions (and polar)']],
            'F':[['A'],['Remove vdW interaction possibility']],
            'I':[['A'],['Remove vdW interaction possibilities']],
            'L':[['A'],['Remove vdW interaction possibilities']],
            'M':[['A'],['Remove vdW interaction possibilities']],
            #'P':[['A'],['Remove vdW interactions / change the shape of the site Normaly Prolines have a structural role and are dangerous to tamper with']],
            'V':[['A'],['Remove vdW interaction possibility']],
            'R':[['A'],['Remove vdW interaction possibility']],
            'N':[['A'],['Remove vdW interaction possibility']],
            'D':[['A'],['Remove vdW interaction possibility']],
            'Q':[['A'],['Remove vdW interaction possibility']],
            'E':[['A'],['Remove vdW interaction possibility']],
            'H':[['A'],['Remove vdW interaction possibility']],
            'K':[['A'],['Remove vdW interaction possibility']],
            'T':[['A'],['Remove vdW interaction possibility']],
            'W':[['A'],['Remove vdW interaction possibility']],
            'Y':[['A'],['Remove vdW interaction possibility']],
            }
        ),
        ('aromatic', { #same for all types of aromatic interactions #fixme check aro_ion
            'Y':[[['L','M'],'A'],['Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities','Remove polar, aromatic and vdW interaction possibilities']],
            'W':[[['L','M'],'A','H'],['Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities','Remove polar, aromatic and vdW interaction possibilities','Weaken aromatic interaction possibilities while retaining polar interaction possibilities']],
            'F':[[['L','M'],'A','Y'],['Remove aromatic interaction possibilities while retaining vdW interaction possibilities','Remove aromatic and vdW interaction possibilities','Prevent Phe-edge to ligand-face aromatic interaction (Should only be attempted when based on a detailed binding mode hypothesis showing the aromatic interaction in question.)']],
            'H':[['A','N',['L','M']],['Remove polar, aromatic and vdW interaction possibilities','Remove aromatic interaction possibilities while retaining polar interaction possibilities','Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities']],
            }
        ),
        ('polar', { #same for all types of polar interactions
            'D':[['L','A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','Remove polar interaction possibility']],
            'E':[[['L','M'],'A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','Remove polar and vdW interaction possibility']],
            'H':[[['L','M'],'A','F'],['Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities','Remove polar, aromatic and vdW interaction possibilities','Remove polar interaction possibility, while retaining aromatic interaction possibilities']],
            'K':[['M','A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','Remove polar and vdW interaction possibility']],
            'N':[['L','A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','Remove polar and vdW interaction possibility']],
            'Q':[[['L','M'],'A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','A - remove the side chain and all interactions']],
            'R':[[['L','M'],'A'],['Remove polar interaction possibility, while retaining vdW interaction possibilities','Remove polar and vdW interaction possibility']],
            'S':[['A'],['Remove polar interaction possibility']],
            'T':[['A','V'],['Remove polar and vdW interaction possibility','Remove polar interaction possibility, while retaining vdW interaction possibilities']],
            'W':[['F','A',['L','M']],['Remove polar interaction possibility, while retaining aromatic interaction possibilities','Remove polar, aromatic and vdW interaction possibilities','Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities']],
            'Y':[['F','A',['L','M']],['Remove polar interaction possibility, while retaining aromatic interaction possibilities (Aromatic interaction cannot be investigated independently from h-bond as none of the polar residues can confromationally match the Tyr OH)','Remove polar, aromatic and vdW interaction possibilities','Remove aromatic and polar interaction possibilities while retaining vdW interaction possibilities']],
            }
        ),
        ('unknown', { #empty
            }
        ),
        ('steric', {
            'A':[[['L','V']],['Increase size to block the binding site']],
            'C':[[['L','I']],['Increase size to block the binding site']],
            'G':[[['A','V','M']],['Add side chain to block the binding site (Glycines are often structurally important and should be mutated with causion - unless closely related receptors show different residues in equivalent position.)']],
            'S':[['L','M'],['Increase size to block the binding site']],
            'T':[['L','M'],['Increase size to block the binding site']],
            'P':[['L','M'],['Increase size to block the binding site (Prolines are often structurally important and should be mutated with causion - unless closely related receptors show different residues in equivalent position.)']],
        })
    ])

AMINO_ACID_GROUP_NAMES = OrderedDict([
        ('hp',     'Hydrophobic - All'),
        ('alhp',   'Hydrophobic - Aliphatic'),
        ('arhp',   'Hydrophobic - Aromatic'),
        ('ar',     'Aromatic'),
        ('pol',    'Polar'),
        ('hb',     'H-Bonding'),
        ('hbd',    'H-Bond Donor'),
        ('hba',    'H-Bond Acceptor'),
        ('charge', 'Charge'),
        ('neg',    'Negative charge'),
        ('pos',    'Positive charge'),
        ('lar',    'Large'),
        ('sma',    'Small'),
        ('any',    'Any feature'),
        ('custom', 'Custom'),
    ])

# Structural rules used to identify different sites
STRUCTURAL_RULES = {    'All': [
                                    # OrderedDict([('Class', 'All'), ('State', 'all'), ('Definition', 'Contract salt-bridge'), ('Generic Position', 'E in upper half of TM, not facing the membrane'), ('Wt AA', 'E'), ('Mut AA', 'D'), ('Note', '')]),
                                    OrderedDict([('Class', 'All'), ('State', 'all'), ('Definition', 'Y7x53 switch removal'), ('Generic Position', '7x53'), ('Wt AA', 'Y'), ('Mut AA', 'L'), ('Note', '')]),
                                    OrderedDict([('Class', 'All'), ('State', 'all'), ('Definition', 'Y7x53 switch removal'), ('Generic Position', '7x53'), ('Wt AA', 'Y'), ('Mut AA', 'A'), ('Note', '')])
                                ],
                        'A': [
                                    OrderedDict([('Class', 'A'), ('State', 'all'), ('Definition', 'Add W3x41'), ('Generic Position', '3x41'), ('Wt AA', 'X'), ('Mut AA', 'W'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '3x49'), ('Wt AA', 'X'), ('Mut AA', 'D'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '3x50'), ('Wt AA', 'X'), ('Mut AA', 'R'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '6x30'), ('Wt AA', 'X'), ('Mut AA', 'E'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) contraction'), ('Generic Position', '3x49'), ('Wt AA', 'E'), ('Mut AA', 'D'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x49'), ('Wt AA', 'D'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x50'), ('Wt AA', 'R'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x50'), ('Wt AA', 'R'), ('Mut AA', 'L'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'all'), ('Definition', 'N6x49 addition'), ('Generic Position', '6x49'), ('Wt AA', 'X'), ('Mut AA', 'N'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '2x50'), ('Wt AA', 'X'), ('Mut AA', 'D'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '3x39'), ('Wt AA', 'X'), ('Mut AA', 'S'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '6x48'), ('Wt AA', 'X'), ('Mut AA', 'W'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '7x49'), ('Wt AA', 'X'), ('Mut AA', 'N'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '2x50'), ('Wt AA', 'X'), ('Mut AA', 'D'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '2x50'), ('Wt AA', 'D/E'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '3x39'), ('Wt AA', 'S/T'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '6x48'), ('Wt AA', 'W'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '7x49'), ('Wt AA', 'N'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '2x50'), ('Wt AA', 'D'), ('Mut AA', 'N'), ('Note', '')]),
                                    OrderedDict([('Class', 'A'), ('State', 'all'), ('Definition', 'Y7x53 addition'), ('Generic Position', '7x53'), ('Wt AA', 'X'), ('Mut AA', 'Y'), ('Note', '')])
                            ],
                        'B': [
                                    OrderedDict([('Class', 'B'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '2x43'), ('Wt AA', 'H'), ('Mut AA', 'A'), ('Note', ' (GPCRdb-B# 2x50)')]),
                                    OrderedDict([('Class', 'B'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x46'), ('Wt AA', 'E'), ('Mut AA', 'A'), ('Note', ' (GPCRdb-B# 3x50)')]),
                                    OrderedDict([('Class', 'B'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '3x39'), ('Wt AA', 'N'), ('Mut AA', 'A'), ('Note', ' (GPCRdb-C# 4x43)')])],
                        'C': [
                                    OrderedDict([('Class', 'C'), ('State', 'all'), ('Definition', 'F7x48 addition'), ('Generic Position', '7x48'), ('Wt AA', 'X'), ('Mut AA', 'F'), ('Note', ' (GPCRdb-B# 7x53)')]),
                                    OrderedDict([('Class', 'C'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '3x50'), ('Wt AA', 'X'), ('Mut AA', 'R'), ('Note', '')]),
                                    OrderedDict([('Class', 'C'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '3x50'), ('Wt AA', 'X'), ('Mut AA', 'K'), ('Note', '')]),
                                    OrderedDict([('Class', 'C'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) addition'), ('Generic Position', '6x33'), ('Wt AA', 'X'), ('Mut AA', 'E'), ('Note', ' (GPCRdb-C# 6x35)')]),
                                    OrderedDict([('Class', 'C'), ('State', 'inactive'), ('Definition', 'Ionic lock (D/ERY) contraction'), ('Generic Position', '6x33'), ('Wt AA', 'E'), ('Mut AA', 'D'), ('Note', ' (GPCRdb-C# 6x35)')]),
                                    OrderedDict([('Class', 'C'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x50'), ('Wt AA', 'R'), ('Mut AA', 'A'), ('Note', '')]),
                                    OrderedDict([('Class', 'C'), ('State', 'active'), ('Definition', 'Ionic lock (D/ERY) removal'), ('Generic Position', '3x50'), ('Wt AA', 'R'), ('Mut AA', 'L'), ('Note', '')]),
                                    OrderedDict([('Class', 'C'), ('State', 'active'), ('Definition', 'Sodium ion site addition'), ('Generic Position', '6x48'), ('Wt AA', 'X'), ('Mut AA', 'W'), ('Note', ' (GPCRdb-C# 6x50)')]),
                                    OrderedDict([('Class', 'C'), ('State', 'inactive'), ('Definition', 'Sodium ion site removal'), ('Generic Position', '6x48'), ('Wt AA', 'W'), ('Mut AA', 'A'), ('Note', ' (GPCRdb-C# 6x50)')]),
                                    OrderedDict([('Class', 'C'), ('State', 'all'), ('Definition', 'Y7x53 addition'), ('Generic Position', '7x53'), ('Wt AA', 'X'), ('Mut AA', 'Y'), ('Note', ' (GPCRdb-B# 7x58)')])
                            ],
                    }

G_PROTEIN_SEGMENTS = OrderedDict([
        ('Full', ['HN','hns1','S1','s1h1','H1','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']),
        ('Structured', ['HN','S1','H1','HA','HB','HC','HD','HE','HF','S2','S3','H2','S4','H3','S5','HG','H4','S6','H5']),
    ])
