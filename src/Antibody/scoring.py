import numpy as np
from itertools import product

from utils.amino_acid import AminoAcids3 as AA3

from Schema.interaction_schema import Interaction, InteractionMember
from typing import List, Literal

def charge_based_heuristic(interactions: List[Interaction]) -> float:

    prop = Literal[
        'negative',
        'positive',
        'polar',
        'hydrophobic',
        'sulfur-bridge',
        'aromatic',
    ]

    order: List[prop] = ['negative', 'positive', 'polar', \
                        'hydrophobic', 'sulfur-bridge', 'aromatic',]

    def typer(m: InteractionMember) -> List[prop]:

        groups = [AA3.CHARGED_NEGATIVE, AA3.CHARGED_POSITIVE, AA3.POLAR, AA3.HYDROPHOBIC, AA3.SULFUR_BRIDGE, AA3.AROMATIC]

        is_in: List[prop] = []

        for group, tag in zip(groups, order):

            if m.residue.object.resname in group:

                is_in.append(tag)

        return is_in

    def scoring(property_a: prop, property_b: prop) -> float:

        loc_a: int = order.index(property_a)
        loc_b: int = order.index(property_b)

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

        r1_types, r2_types = typer(r1), typer(r2)

        combs = product(r1_types, r2_types)

        for c in combs:

            total += scoring(*c)

    return total