from pydantic import BaseModel

from typing import List, Optional

from Schema.structure_schema import Protien, Residue, Atom

class InteractionMember(BaseModel):

    protien: Protien
    residue: Residue
    atom: Optional[Atom]

class Interaction(BaseModel):

    members: List[InteractionMember]

class CentralCarbonDistance(Interaction):

    central_carbon_distance: float

class AtomClosestDistance(CentralCarbonDistance):

    closest_atom_distance: float
