from pydantic import BaseModel, ConfigDict, model_validator, FilePath

from Bio.PDB.Structure import Structure as BioStructure
from Bio.PDB.Chain import Chain as BioChain
from Bio.PDB.Residue import Residue as BioResidue
from Bio.PDB.Atom import Atom as BioAtom, DisorderedAtom

# from Schema import Interaction, Causes circular import

from typing import Tuple, Dict, List, Any, Literal

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

    # Disorderd Atoms act the same as Atoms
    object: BioAtom | DisorderedAtom

    id: AtomFullID

    model_config = ConfigDict(arbitrary_types_allowed=True)

class Residue(BaseModel):

    # points to bio object
    object: BioResidue

    id: ResidueFullID

    atoms: List[Atom]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def ca(self) -> Atom:
        """Gets central carbon atom"""

        for a in self.atoms:

            if a.object.id == 'CA':

                return a
        
        # Sometimes do not have a CA, ex ABU
        return self.atoms[0]

class Protien(BaseModel):

    # points to Bio object
    object: BioChain

    # Does not have to be ordered
    residues: List[Residue]

    model_config = ConfigDict(arbitrary_types_allowed=True)

class Structure(BaseModel):

    object: BioStructure
    file: FilePath
    file_type: Literal['pdb', 'cif']

    participants: List[Protien]
    interactions: List[Any] # List of Interaction

    model_config = ConfigDict(arbitrary_types_allowed=True)