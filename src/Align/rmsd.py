from Bio.PDB import MMCIFParser, Superimposer
import pandas as pd

def rmsd(file_a: str, file_b: str):
    parser = MMCIFParser(QUIET=True)
    ref = parser.get_structure("ref", file_a)
    mob = parser.get_structure("mob", file_b)

    # Extract Cα atoms
    ref_atoms = [a for a in ref.get_atoms() if a.name == "CA"]
    mob_atoms = [a for a in mob.get_atoms() if a.name == "CA"]

    sup = Superimposer()
    sup.set_atoms(ref_atoms[:len(mob_atoms)], mob_atoms[:len(ref_atoms)])
    sup.apply(mob.get_atoms())

    print("RMSD:", sup.rms)
    return sup.rms

if __name__ == '__main__':

    df = pd.DataFrame(columns = ['s1', 's2', 'rmsd'])

    structs = ['reelin56_apoer', *sum([[s + str(i) for i in range(5)] for s in ['no_ca_', 'ca_']], [])]

    for i, s1 in enumerate(structs):
        for s2 in structs[i:]:
            print(s1, s2)
            r = rmsd(s1 + '.cif', s2 + '.cif')

            df.loc[len(df)] = [s1, s2, str(r)]

    df.to_csv('af3_rmsd_results.csv')