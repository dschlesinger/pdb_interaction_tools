from .distance import residue_interaction_by_ca, residue_closest_atomic_interaction, residue_by_id
from .AminoAcids import AminoAcids3
from .ViewTools import ViewHandler

### Create pseudonyms for function

# residue_interaction_by_distance
ribc = by_residue_ca =  residue_interaction_by_ca
rcai = by_residue_closest = residue_closest_atomic_interaction