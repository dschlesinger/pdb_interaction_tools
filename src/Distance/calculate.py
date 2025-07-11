import numpy as np

from Schema import Residue, Protien, Interaction, InteractionMember, CentralCarbonDistance, AtomClosestDistance

from typing import List

class CalculateInteraction:

    def interactions_by_member(r: Residue, interactions: List[Interaction]) -> List[Interaction]:

        found = []

        for i in interactions:

            a, b, *_ = i.members

            if r in [a.residue, b.residue]:

                found.append(i)

        return found

    @staticmethod
    def by_closest_atoms(A: Protien, B: Protien, central_atom_cutoff: float = 10.0, atomic_cutoff: float = 6.0) -> List[AtomClosestDistance]:

        filtered: List[CentralCarbonDistance] = \
            CalculateInteraction.by_central_carbon(A, B, central_atom_cutoff)

        atomic_interactions: List[AtomClosestDistance] = []

        for interaction in filtered:

            # Ignore if more than 2
            r1, r2, *_ = interaction.members

            r1_coords = np.array([atom.object.coord for atom in r1.residue.atoms])
            r2_coords = np.array([atom.object.coord for atom in r2.residue.atoms])

            # Rows are r1 residues, Columns are B residues
            diff = r1_coords[:, np.newaxis, :] - r2_coords[np.newaxis, :, :]

            # Uses euclidean distance
            distance_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

            if distance_matrix.min() < atomic_cutoff:

                r1_pos, r2_pos = np.unravel_index(distance_matrix.argmin(), distance_matrix.shape)

                r1.atom = r1.residue.atoms[r1_pos]
                r2.atom = r2.residue.atoms[r2_pos]
                atomic_interactions.append(
                    AtomClosestDistance(
                        members=interaction.members,
                        central_carbon_distance=interaction.central_carbon_distance,
                        closest_atom_distance=distance_matrix.min().item(),
                    )
                )

        return atomic_interactions


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