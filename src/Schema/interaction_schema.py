from pydantic import BaseModel, model_validator, Field, ConfigDict

from Bio.PDB.Structure import Structure as BioStructure
from Bio.PDB.Chain import Chain as BioChain
from Bio.PDB.Residue import Residue as BioResidue
from Bio.PDB.Atom import Atom as BioAtom, DisorderedAtom

from typing import List, Tuple, Dict, Annotated

class ResidueFullID(BaseModel):

    full_id: Tuple

    # Basiclly PDB code
    structure_id: str

    # Used in structure[int]
    model_id: int

    # ex 'A' or 'B'
    chain_id: str

    # Residue ID tuple (' ', int, ' ')
    hetero_field: str # says water or smth
    sequence_id: int
    intersection_id: str # if is from other chain, should not be a problem

    # As model is passed in as tuple
    @model_validator(mode='before')
    @classmethod
    def process_full_id(cls, data: Dict[str, str | Tuple]) -> Dict[str, int | str]:
        """example {full_id: ('9fg1', 0, 'B', (' ', 18, ' ')), list_position: int}"""

        kw = data

        kw['structure_id'], \
            kw['model_id'], \
            kw['chain_id'], \
            rid = data['full_id']
        
        kw['hetero_field'], \
        kw['sequence_id'], \
        kw['intersection_id'] = rid

        return kw

class AtomFullID(ResidueFullID):

    alternative_location_id: str # only if mutliple possible locations

class Atom(BaseModel):

    object: BioAtom | DisorderedAtom

    id: AtomFullID

    model_config = ConfigDict(arbitrary_types_allowed=True)

class Residue(BaseModel):

    # points to bio object
    object: BioResidue

    id: ResidueFullID

    atoms: List[Atom]

    model_config = ConfigDict(arbitrary_types_allowed=True)

class Protien(BaseModel):

    # points to Bio object
    object: BioChain

    # Does not have to be ordered
    amino_acids: List[Residue]

    model_config = ConfigDict(arbitrary_types_allowed=True)

class InteractionMember(BaseModel):

    protien: Protien
    residue: Residue

class Interaction(BaseModel):

    members: List[InteractionMember]

class CentralCarbonDistance(Interaction):

    central_carbon_distance: float

class AtomClosestDistance(CentralCarbonDistance):

    closest_atom_distance: float


class Structure(BaseModel):

    object: BioStructure

    participants: List[Protien]
    interactions: List[Interaction]

    model_config = ConfigDict(arbitrary_types_allowed=True)
