from ccdc import search
import pandas as pd
from io import StringIO


def smarts_to_hits(
    smarts,
    db,
    atom_property=False,
    max_hits_per_structure=None,
    torsional_constraint=None,
):
    searcher = search.SubstructureSearch()
    searcher.Settings(max_hits_per_structure=max_hits_per_structure)
    smarts = search.SMARTSSubstructure(smarts)
    if atom_property:
        smarts.atoms[0].label_match = "^_Z"
    searcher.add_substructure(smarts)

    if torsional_constraint:
        searcher.add_torsion_angle_constraint(*torsional_constraint)
    return searcher.search(db, max_hits_per_structure=max_hits_per_structure)


def ccdc_mol_to_gauss_df(ccdc_mol):
    atom_labels = []
    xs = []
    ys = []
    zs = []

    for a in ccdc_mol.atoms:
        atom_labels.append(a.label)
        xs.append(a.coordinates[0])
        ys.append(a.coordinates[1])
        zs.append(a.coordinates[2])
    df = pd.DataFrame({"atom_label": atom_labels, "x": xs, "y": ys, "z": zs})
    return df


def mol2_to_df(mol2_string):
    record_types = mol2_string.split("@<TRIPOS>")
    atoms_string = "\n".join(record_types[2].split("\n")[1:])
    bonds_string = "\n".join(record_types[3].split("\n")[1:])
    mol_string = "\n".join(record_types[1].split("\n")[1:])
    for atom_label in [
        "_ZN",
        "_UN",
        "_ZO",
        "_UO",
        "_ZC",
        "_UC",
        "_ZCl",
        "_UCl",
    ]:

        atoms_string = atoms_string.replace(f"{atom_label} ", atom_label)
    atoms_df = pd.read_csv(StringIO(atoms_string), sep="\s+", header=None)

    return atoms_df, bonds_string, mol_string
