from utils.load_pdb import load_structure
from Distance.calculate import CalculateInteraction
from Antibody.scoring import charge_based_heuristic

def test_atomic_distance() -> None:

    s = load_structure('9fg1.pdb')

    a, b, *_ = s.participants

    c = CalculateInteraction.by_closest_atoms(a, b)

    for ci in c:

        a, b, *_ = ci.members

        s = charge_based_heuristic([ci])

        print(f'{a.residue.object.resname} -> {b.residue.object.resname}: {s}')
