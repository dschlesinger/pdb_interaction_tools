# Just import it
from utils.visualize import show_antibody_antigen
from utils.load_pdb import load_structure
from Distance.calculate import CalculateInteraction

from utils.amino_acid import AminoAcids3 as AA3

def test_run() -> None: 
    s = load_structure('9fg1.pdb')

    x, y, *_ = s.participants

    i = CalculateInteraction.by_closest_atoms(x, y)

    a = i[0].members[0]

    print(AA3.color_by_charge(a), a.residue)

    show_antibody_antigen(s, highlight=i)

    