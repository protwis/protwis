from collections import OrderedDict

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
        ('J', 'Xle'),
        ('X', 'Xaa'),
    ])

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
            'A':[[['L','V','I']],['Increase size to "block" binding']],
            'C':[[['L','V','I']],['Increase size to "block" binding']],
            'F':[['A'],['Remove vdW interactions / change the shape of the site']],
            'I':[['A'],['Remove vdW interactions / change the shape of the site']],
            'L':[['A'],['Remove vdW interactions / change the shape of the site']],
            'M':[['A'],['Remove vdW interactions / change the shape of the site']],
            'P':[['A'],['Remove vdW interactions / change the shape of the site Normaly Prolines have a structural role and are dangerous to tamper with']],
            'V':[['A'],['Remove vdW interactions / change the shape of the site']],
            'W':[['A'],['Remove vdW interactions / change the shape of the site']],
            'Y':[['A'],['Remove vdW interactions / change the shape of the site']]
            }
        ),
        ('aromatic', { #same for all types of aromatic interactions #fixme check aro_ion 
            'Y':[[['L','M'],'A'],['L,M - Remove aromaticity - keep the size','A - remove the side chain and all interactions']],
            'W':[[['L','M'],'A'],['L,M - Remove aromaticity - keep the size','A - remove the side chain and all interactions']],
            'F':[[['L','M'],'A','Y'],['L,M - Remove aromaticity - keep the size','A - remove the side chain and all interactions','Y - OH may prevent the aromatic interaction']],
            'H':[[['L','M'],'A','F'],['L,M - Remove aromaticity - keep the size','A - remove the side chain and all interactions','F - Make it more aromatic (stronger interaction)']],
            }
        ),
        ('polar', { #same for all types of polar interactions
            'D':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'E':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'H':[[['L','M'],'A','F'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions','F - remove h-bond and keep aromatic']],
            'K':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'N':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'Q':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'R':[[['L','M'],'A'],['L,M - Remove interaction - keep some size - M more flexible','A - remove the side chain and all interactions']],
            'S':[['A','L'],['A - remove the side chain and all interactions','L - introduce bulk and block the site']],
            'T':[['A','V','L'],['A - remove the side chain and all interactions','V - remove the h-bond and keep the size','L - introduce bulk and block the site']],
            'W':[['F','A'],['F - remove h-bond and keep aromaticity and size','A - remove the side chain and all interactions']],
            'Y':[['F','A'],['F - remove h-bond and keep aromaticity and size','A - remove the side chain and all interactions']],
            }
        ),
        ('unknown', { #empty
            }
        )
    ])