from Bio.PDB import PDBParser, MMCIFParser, Structure as BioStructure, PDBExceptions
from Schema import Atom, Residue, Protien, Structure, AtomFullID, ResidueFullID

from pydantic import FilePath
from typing import Literal, Optional
import os

import warnings


def load_structure(file_path: FilePath, file_type: Optional[Literal['pdb', 'cif']] = None, frame: Optional[int] = 0) -> Structure:
    """Loads file and returns structure, takes pdb or cif
    Sets structure id to the filename before the first .

    Args:
        file_path (FilePath): path to file, if no file_type infers from 9fg1.pdb <--
        file_type (Optional[Literal[&#39;pdb&#39;, &#39;cif&#39;]]): pdb | cif
        frame (Optional[int]): what frame of structure to use for Protiens, default is 0

    Raises:
        FileNotFoundError: Could not find given file
        ValueError: if file_type not given and filename.___ is not pdb or cif

    Returns:
        Structure: structure as Schema.interaction_schema.Structure
    """    

    if not os.path.exists(file_path):

        raise FileNotFoundError(file_path)

    file_id, *_, file_ending = file_path.split('.')

    file_type = file_type or file_ending

    match file_type:

        case 'pdb':

            parser = PDBParser()

        case 'cif':

            parser = MMCIFParser()
        
        case _:

            raise ValueError(f'{file_type} is not pdb or cif')
    
    # So doesn't effect global settings
    with warnings.catch_warnings():

        # Gets annoying sometimes
        warnings.filterwarnings("ignore", category=PDBExceptions.PDBConstructionWarning)

        bio_struct = parser.get_structure(file_id, file_path)

    structure = Structure(
        participants=[],
        interactions=[],
        object=bio_struct,
    )

    # Iterate over subparts
    for chain in structure.object[frame]:

        c = Protien(
            object=chain,
            residues=[],
        )

        structure.participants.append(c)

        for residue in c.object:

            r = Residue(
                object=residue,
                id=ResidueFullID(
                    full_id=residue.full_id
                ),
                atoms=[],
            )  

            c.residues.append(r)

            for atom in r.object:

                a = Atom(
                    object=atom,
                    id=AtomFullID(
                        full_id=atom.full_id[:-1],
                        alternative_location_id=atom.full_id[-1][-1]
                    )
                )

                r.atoms.append(a)
    
    return structure

    