import numpy as np
from typing import List

"""
Data source v

@unknown{unknown,
    author = {Zhang, Grace and Su, Zhaoqian and Zhang, Tom and Wu, Yinghao},
    year = {2023},
    month = {12},
    pages = {},
    title = {Machine-learning-based Structural Analysis of Interactions between Antibodies and Antigens},
    journal = {bioRxiv : the preprint server for biology},
    doi = {10.1101/2023.12.06.570397}
}
"""

class AAPairFrequency:
    """Wrapper class for AB AG AA pair frequency at binding interface

    AA_ORDER: List[str] - 3 letter AA's, index in this list is the row and column in the matrix
    
    """    

    AA_ORDER: List[str] = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    FREQUENCY_MATRIX: np.ndarray = np.load('residue_frequency.py')