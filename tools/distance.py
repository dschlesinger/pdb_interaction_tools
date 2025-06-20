import numpy as np
from Bio import PDB

from typing import List, Dict, Tuple, Self
from dataclasses import dataclass

import json

@dataclass
class Interaction:

  id: str
  distance: float

@dataclass
class InteractingResidue:

  chain: str
  position: int # in chain

  # to keep track of interactions
  id: str

  interactions: List[Interaction]

  @classmethod
  def from_json(cls, instance: Dict) -> Self:
      
      interactions = [Interaction(**i) for i in instance['interactions']]

      del instance['interactions']
      
      return cls(**instance, interactions=interactions)

  def to_json(self) -> str:

    base_json = self.__dict__

    # __dict__ interactions
    base_json['interactions'] = [i.__dict__ for i in base_json['interactions']]

    return base_json

def residue_interaction_by_distance(a: PDB.Chain.Chain, b: PDB.Chain.Chain, distance: float = 6.0) -> List[InteractingResidue]:
    """Returns a list of all residues with at least one interaction under the distance limit

    Args:
        a (PDB.Chain.Chain): First protien
        b (PDB.Chain.Chain): Second protien
        distance (float): in angstroms

    Returns:
        List[InteractingResidue]: List of all residues with a hit
    """  

    # Get central carbon atoms
    a_ca = [atom for atom in a.get_atoms() if atom.id == 'CA']
    b_ca = [atom for atom in b.get_atoms() if atom.id == 'CA']

    # Get coords of these carbon atoms
    a_mat = np.array([atom.coord for atom in a_ca], dtype=np.float32)
    b_mat = np.array([atom.coord for atom in b_ca], dtype=np.float32)

    # Rows are A residues, Columns are B residues
    diff = a_mat[:, np.newaxis, :] - b_mat[np.newaxis, :, :]

    # Uses euclidean distance
    distance_mat = np.sqrt(np.sum(diff ** 2, axis=2))

    # Get hits
    mask = distance_mat < distance

    below_cutoff = distance_mat[mask]

    # key is f"{chain.id}:{position}"
    interacting_residues: Dict[str, InteractingResidue] = {}

    for idx, euc_dist in zip(np.argwhere(mask), below_cutoff):

      # Correct idx to residue position in parent
      pos = (a_ca[idx[0]].get_parent().get_id()[1], b_ca[idx[1]].get_parent().get_id()[1])

      # Keys for dictionary
      a_key = ':'.join((a.id, str(pos[0])))
      b_key = ':'.join((b.id, str(pos[1])))

      # Create Interaction objects
      a_interaction = Interaction(id=a_key, distance=float(euc_dist))
      b_interaction = Interaction(id=b_key, distance=float(euc_dist))

      # Add to chain a
      if a_key in interacting_residues:
            interacting_residues[a_key].interactions.append(b_interaction)
      else:
          interacting_residues[a_key] = InteractingResidue(
              chain=a.id,
              position=pos[0],
              id=a_key,
              interactions=[b_interaction]
          )

      # Add to chain b
      if b_key in interacting_residues:
            interacting_residues[b_key].interactions.append(a_interaction)
      else:
          interacting_residues[b_key] = InteractingResidue(
              chain=b.id,
              position=pos[1],
              id=b_key,
              interactions=[a_interaction]
          )

    return list(interacting_residues.values())