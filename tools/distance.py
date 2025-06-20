import numpy as np
from Bio import PDB

from typing import List, Dict, Tuple, Any
from dataclasses import dataclass

import json

@dataclass
class AtomicID:
    """Identifier for molecules in the protien"""
    chain: str

    residue_number: int
    residue_type: str

    atom_type: str
    atom_number: int

@dataclass
class Interaction:

  id: str
  global_distance: float

  # If using by closest atomic interaction
  closest_atomic_diffence: float | None = None
  self_closest: AtomicID | None = None
  other_closest: AtomicID | None = None

@dataclass
class InteractingResidue:

  chain: str
  position: int # in chain

  # For later analysis, in 3 letter format
  residue: str

  # to keep track of interactions
  id: str

  interactions: List[Interaction]

  @classmethod
  def from_json(cls, instance: Dict):

      for interaction in instance['interactions']:
          
        if interaction['self_closest']:

            interaction['self_closest'] = AtomicID(**interaction['self_closest'])
        
        if interaction['other_closest']:

            interaction['other_closest'] = AtomicID(**interaction['other_closest'])
      
      interactions = [Interaction(**i) for i in instance['interactions']]

      del instance['interactions']
      
      return cls(**instance, interactions=interactions)

  def to_json(self) -> str:

    base_json = self.__dict__

    # __dict__ interactions
    for interaction in base_json['interactions']:

        if interaction.self_closest:

            interaction.self_closest = interaction.self_closest.__dict__

        if interaction.other_closest:

            interaction.other_closest = interaction.other_closest.__dict__

    base_json['interactions'] = [i.__dict__ for i in base_json['interactions']]

    return base_json

def residue_by_id(interacting_residue: List[InteractingResidue], id: str) -> InteractingResidue | None:
    """Returns None if not found"""

    for ir in interacting_residue:

        if ir.id == id:

            return ir

def Bio_object_by_residue(a: PDB.Chain.Chain, r: InteractingResidue) -> Any: # idk what type

    # Get residue
    for res in a.get_residues():

        if r.position == res.id[1]:

            return res

def residue_closest_atomic_interaction(a: PDB.Chain.Chain, b: PDB.Chain.Chain, distance: float = 6.0) -> List[InteractingResidue]:

    by_CA_distance: List[InteractingResidue] = residue_interaction_by_ca(a, b, distance)

    # Fill in None fields
    for self_residue in by_CA_distance:

        if self_residue.chain == a.id:
            self_chain = a
            other_chain = b

        else:
            self_chain = b
            other_chain = a

        self_Bio_object = Bio_object_by_residue(self_chain, self_residue)

        self_atom_coord_matrix = np.array([atom.coord for atom in self_Bio_object.get_atoms()])

        for interaction in self_residue.interactions:

            other_residue: InteractingResidue = residue_by_id(by_CA_distance, interaction.id)

            other_Bio_object = Bio_object_by_residue(other_chain, other_residue)

            other_atom_coord_matrix = np.array([atom.coord for atom in other_Bio_object.get_atoms()])

            # Self is rows other is columns
            diff = self_atom_coord_matrix[:, np.newaxis, :] - other_atom_coord_matrix[np.newaxis, :, :]

            # Uses euclidean distance
            distance_mat = np.sqrt(np.sum(diff ** 2, axis=2))

            closest_atoms = distance_mat.argmin()

            # Get coords from flat index
            closest_atoms = np.unravel_index(closest_atoms, distance_mat.shape)

            interaction.closest_atomic_diffence = distance_mat[closest_atoms].item()

            self_atom = list(self_Bio_object.get_atoms())[closest_atoms[0].item()]
            other_atom = list(other_Bio_object.get_atoms())[closest_atoms[1].item()]

            interaction.self_closest = AtomicID(
                chain=self_chain.id,
                residue_number=self_residue.position,
                residue_type=self_residue.residue,
                atom_number=closest_atoms[0].item(),
                atom_type=self_atom.element,
            )

            interaction.other_closest = AtomicID(
                chain=other_chain.id,
                residue_number=other_residue.position,
                residue_type=other_residue.residue,
                atom_number=closest_atoms[1].item(),
                atom_type=other_atom.element,
            )
    
    # bc why not
    by_atomic_diffrence = by_CA_distance

    return by_atomic_diffrence


def residue_interaction_by_ca(a: PDB.Chain.Chain, b: PDB.Chain.Chain, distance: float = 6.0) -> List[InteractingResidue]:
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
      a_residue = a_ca[idx[0]].get_parent()
      b_residue = b_ca[idx[1]].get_parent()

      pos = (a_residue.get_id()[1], b_residue.get_id()[1])

      # Keys for dictionary
      a_key = ':'.join((a.id, str(pos[0])))
      b_key = ':'.join((b.id, str(pos[1])))

      # Create Interaction objects
      a_interaction = Interaction(id=a_key, global_distance=float(euc_dist))
      b_interaction = Interaction(id=b_key, global_distance=float(euc_dist))

      # Add to chain a
      if a_key in interacting_residues:
            interacting_residues[a_key].interactions.append(b_interaction)
      else:
          interacting_residues[a_key] = InteractingResidue(
              chain=a.id,
              position=pos[0],
              residue=a_residue.get_resname(),
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
              residue=b_residue.get_resname(),
              id=b_key,
              interactions=[a_interaction]
          )

    return list(interacting_residues.values())