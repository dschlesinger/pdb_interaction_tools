class AminoAcids3:
    """Class stores 3 letter amino acid residue names and groups"""

    # TRP is weird because it has N-H but is large so kinda non-polar
    # So I put in both

    # Amino Acid groups
    ALL = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL', 'SEC'  # Including Selenocysteine
    ]

    HYDROPHOBIC = [
        'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'MET', 'PHE', 'TRP'
    ]

    CHARGED = [
        'ARG', 'LYS', 'HIS',  # Positively charged at physiological pH
        'ASP', 'GLU'          # Negatively charged
    ]

    CHARGED_POSITIVE = [
        'ARG', 'LYS', 'HIS'
    ]

    CHARGED_NEGATIVE = [
        'ASP', 'GLU'
    ]

    POLAR = [
        'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS', 'HIS', 'TRP', 'SEC'
    ]

    SULFUR_BRIDGE = [
        'CYS', 'SEC'  # Selenocysteine is the selenium analog of cysteine
    ]

    AROMATIC = [
        'PHE', 'TYR', 'TRP', 'HIS'
    ]
    # ChatGPT Generated check with wikipedia