import numpy as np
from itertools import product

from utils.amino_acid import AminoAcids3 as AA3

from Schema import Interaction, InteractionMember
from typing import List, Literal

def charge_based_heuristic(interactions: List[Interaction]) -> float:

    def scoring(property_a: AA3.charge_property, property_b: AA3.charge_property) -> float:

        loc_a: int = AA3.order.index(property_a)
        loc_b: int = AA3.order.index(property_b)

        # Upper right triangle
        score_matrix: np.ndarray = np.array(
        [
            [-2, 2, 1, -1, 0, 0], # Negative
            [0, -2, 1, -1, 0, 0], # Positive
            [0, 0, 1, -1, 0, 0], # Polar
            [0, 0, 0, 1, 0, 1], # Hydrophobic
            [0, 0, 0, 0, 4, 0], # Sulfur bridge
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

            total += scoring(*c)

    return total