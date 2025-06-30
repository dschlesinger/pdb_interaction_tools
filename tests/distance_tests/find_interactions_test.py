from utils.load_pdb import load_structure
from Distance.calculate import CalculateInteraction
from Antibody.scoring import charge_based_heuristic

def test_central_carbon() -> None:

    s = load_structure('9fg1.pdb')

    a, b, *_ = s.participants

    c = CalculateInteraction.by_central_carbon(a, b)

    # print(c[0].central_carbon_distance)
    # print(len(c))

def test_atomic_distance() -> None:

    s = load_structure('9fg1.pdb')

    a, b, *_ = s.participants

    c = CalculateInteraction.by_closest_atoms(a, b)

    print([ci.atom for ci in c[0].members])

def test_charged_based_heuristic() -> None:

    s = load_structure('9fg1.pdb')

    a, b, *_ = s.participants

    c = CalculateInteraction.by_closest_atoms(a, b)

    print(charge_based_heuristic(c))