from Schema import ResidueFullID, AtomFullID

def test_residue_full_id() -> None:

    full_id = ('9fg1', 0, 'B', (' ', 0, ' '))
    list_position = 0

    ResidueFullID(full_id=full_id, list_position=list_position)

def test_atom_full_id() -> None:

    full_id = ('9fg1', 0, 'B', (' ', 0, ' '), ('N', ' '))
    list_position = 0

    _, altloc = full_id[-1]

    full_id = full_id[:-1]

    AtomFullID(full_id=full_id, list_position=list_position, alternative_location_id=altloc)