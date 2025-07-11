import py3Dmol

from Schema import Interaction, Structure, InteractionMember

from .amino_acid import AminoAcids3 as AA3

from typing import List, Literal, List

def show_antibody_antigen(
    s: Structure,
    *,
    highlight: List[Interaction] = [],
    highlight_by: Literal['charge'] = 'charge',
    default_color: str = 'black',
) -> py3Dmol.view:
    
    with open(s.file, 'r') as f:

        protien_stuff = f.read()

    view = py3Dmol.view(width=800, height=800)

    view.addModel(protien_stuff, s.file_type)

    view.setStyle({'cartoon': {'color': default_color}})
    
    residues_to_color: List[InteractionMember] = sum(
            [[im for im in interaction.members] for interaction in highlight]
            , []
        )
    
    for r in residues_to_color:

        color = AA3.color_by_charge(r, default_color=default_color)

        view.setStyle({'chain': r.protien.object.id, 'resn': r.residue.object.get_id()[1]}, {'stick': {'color': color}})

    return view