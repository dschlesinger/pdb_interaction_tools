from utils.load_pdb import load_structure

from Schema.structure_schema import Mutation

s = None

def test_structure_loading() -> None:
    global s
    s = load_structure('9fg1.pdb')

def test_single_letter() -> None:

    a, b, *_ = s.participants

    r = a.residues[0]

    print(r.object.resname, r.single_letter_resname)

def test_sequence() -> None:

    a, b, *_ = s.participants

    m = Mutation(
        chain=a.object.id,
        resi=2,
        to='A',
    )

    print(a.seq()[:5], a.seq([m])[:5])

