# Functions to score AB AG binding and specific AA pairs by contribution

import numpy as np

if __name__ == '__main__':
    import sys
    sys.path.append('/workspaces/pdb_interaction_tools')

from tools import AminoAcids3 as AA3
from tools.distance import InteractingResidue, residue_by_id

# Relitave frequency of AA pairs in AB AG interfaces
# from tools.antibody.residue_frequency import AAPairFrequency

from typing import List, Set, Tuple, Dict

def property_based_heuristic(interactions: List[InteractingResidue]) -> int:

    def scoring(a_prop: List[str], b_prop: List[str]) -> int:
        """Scores individual pairs"""

        # Alphabetical for the tuple
        SCORE_MAP: Dict[Tuple[str, str], int] = {
            ('aromatic', 'aromatic'): 4,
            ('aromatic', 'negative'): 0,
            ('aromatic', 'positive'): 0,
            ('aromatic', 'polar'): 0,
            ('aromatic', 'hydrophobic'): 1,
            ('hydrophobic', 'negative'): -1,
            ('hydrophobic', 'positive'): -1,
            ('hydrophobic', 'polar'): -1,
            ('hydrophobic', 'hydrophobic'): 1,
            ('negative', 'negative'): -2,
            ('negative', 'positive'): 2,
            ('negative', 'polar'): 1,
            ('polar', 'polar'): 1,
            ('polar', 'positive'): 1,
            ('positive', 'positive'): -2,
        }

        total: int = 0

        for a in a_prop:
            for b in b_prop:

                total += SCORE_MAP[(min(a, b), max(a, b))]
        
        return total

    # storers hashes of interactions evaluated
    evaluated_interactions: Set[int] = set()
    score: int = 0

    TAG_LIST_MAP = [
        ('hydrophobic', AA3.HYDROPHOBIC),
        ('polar', AA3.POLAR),
        ('negative', AA3.CHARGED_NEGATIVE),
        ('positive', AA3.CHARGED_POSITIVE),
        ('aromatic', AA3.AROMATIC),
    ]

    for residue in interactions:

        for interaction in residue.interactions:

            if interaction.hash in evaluated_interactions:

                continue

            evaluated_interactions.add(interaction.hash)

            a_prop = []
            b_prop = []

            for tag, l in TAG_LIST_MAP:

                if residue.residue in l:
                    a_prop.append(tag)

                if (r := residue_by_id(interactions, interaction.id)) and r.residue in l:            
                    b_prop.append(tag)

            score += scoring(a_prop, b_prop)

    return score


if __name__ == '__main__':

    import Bio
    from Bio.PDB import PDBParser
    from Bio.PDB.Chain import Chain

    from tools import residue_closest_atomic_interaction as rcai

    import warnings

    warnings.filterwarnings("ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning)

    parser = PDBParser()
    structure = parser.get_structure("9fg1", "/workspaces/pdb_interaction_tools/9fg1.pdb")

    A: Chain = structure[0]['A']
    B: Chain = structure[0]['B']

    interactions = rcai(A, B, distance=10.0, atomic_distance=6.0)

    print(property_based_heuristic(interactions))