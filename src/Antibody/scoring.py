import numpy as np
from itertools import product

from utils.amino_acid import AminoAcids3 as AA3

from Schema import Interaction, AtomClosestDistance, InteractionMember
from typing import List, Literal

def charge_based_heuristic(interactions: List[AtomClosestDistance], cutoff: float = 6) -> float:

    def scoring(property_a: AA3.charge_property, property_b: AA3.charge_property) -> float:

        # Order ['hydrophobic', 'polar', 'negative', \
        #            'positive', 'sulfur-bridge', 'aromatic',]
        loc_a: int = AA3.order.index(property_a)
        loc_b: int = AA3.order.index(property_b)

        # Upper right triangle
        score_matrix: np.ndarray = np.array(
        [
            [1, -1, -1, -1, 0, 1], # Hydrophobic
            [0, 1, 1, 1, 0, 0], # Polar
            [0, 0, -2, 2, 0, 0], # Negative
            [0, 0, 0,-2, 0, 0], # Positive
            [0, 0, 0, 0, 4, 0], # Sulfur Bridge
            [0, 0, 0, 0, 0, 4], # Aromatic
        ],
            dtype=np.float32,
        )

        return score_matrix[min(loc_a, loc_b), max(loc_a, loc_b)]
    
    total: float = 0

    for interaction in interactions:

        r1, r2, *_ = interaction.members

        r1_types, r2_types = AA3.charge_typer(r1), AA3.charge_typer(r2)

        combs = product(r1_types, r2_types)

        for c in combs:

            # Score * distance
            s = scoring(*c) * (1 - ((interaction.closest_atom_distance / cutoff) ** 2))

            total += s

    return total