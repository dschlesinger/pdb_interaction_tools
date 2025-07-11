from Schema import InteractionMember

from typing import List, Literal, Dict

CHARGE_COLOR_MAP: Dict[str, str] = {
    'aromatic': 'green',
    'negative': 'blue',
    'positive': 'red',
    'polar': 'white',
    'hydrophobic': 'grey',
    'sulfur-bridge': 'yellow',
}

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

    charge_property = Literal['polar', 'hydrophobic', 'sulfur-bridge', \
                              'negative', 'positive', 'aromatic',]

    order: List[charge_property] = ['polar', 'hydrophobic', 'sulfur-bridge', \
                              'negative', 'positive', 'aromatic',]

    @classmethod
    def color_by_charge(cls, m: InteractionMember, default_color: str = 'black') -> str:

        # Default and Glyicne
        color = default_color

        properties = cls.charge_typer(m)

        for p in properties:

            color = CHARGE_COLOR_MAP[p]

        return color


    @classmethod
    def charge_typer(cls, m: InteractionMember) -> List['charge_property']:

        groups = [cls.POLAR, cls.HYDROPHOBIC, cls.CHARGED_NEGATIVE, cls.CHARGED_POSITIVE, cls.SULFUR_BRIDGE, cls.AROMATIC]

        is_in: List['charge_property'] = []

        for group, tag in zip(groups, cls.order):

            if m.residue.object.resname in group:

                is_in.append(tag)

        return is_in