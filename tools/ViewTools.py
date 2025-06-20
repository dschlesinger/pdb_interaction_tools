from typing import List, Any

class ViewHandler:

    @staticmethod
    def highlight_aa(view: Any, aa_code: List[str], color: str) -> None:
        """Takes py3dmol view and highlights any amino acids that match with the given color"""

        for aa in aa_code:

            view.setStyle({'resn': aa}, {'cartoon': {'color': color}})