import numpy as np

from Schema.interaction_schema import Protien, Interaction, InteractionMember, CentralCarbonDistance, AtomClosestDistance

from typing import List

class CalculateInteraction:

    @staticmethod
    def by_central_carbon(A: Protien, B: Protien, distance_cutoff: float = 6.0) -> List[CentralCarbonDistance]:

        A_coords = np.array([r.ca().object.coord for r in A.residues])
        B_coords = np.array([r.ca().object.coord for r in B.residues])

        # Rows are A residues, Columns are B residues
        diff = A_coords[:, np.newaxis, :] - B_coords[np.newaxis, :, :]

        # Uses euclidean distance
        distance_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

        # Mask and argwhere
        mask = distance_matrix < distance_cutoff
        below_cutoff = distance_matrix[mask]

        interactions: List[CentralCarbonDistance] = []

        for idx, distance in zip(np.argwhere(mask), below_cutoff):

            A_residue = InteractionMember(
                protien=A,
                residue=A.residues[idx[0]],
                atom=None,
            )

            B_residue = InteractionMember(
                protien=B,
                residue=B.residues[idx[1]],
                atom=None,
            )

            i: Interaction = CentralCarbonDistance(
                members=[A_residue, B_residue],
                central_carbon_distance=distance.item(),
            )

            interactions.append(i)

        return interactions