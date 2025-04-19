import argparse
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Generate conformers with fixed fragment.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument("-l", "--ligand", help="Ligand file.", default="ligand.smi")

    parser.add_argument(
        "-lf",
        "--ligand_format",
        help="Ligand SMILES.",
        choices=["sdf", "smi"],
        default="smi",
    )

    parser.add_argument(
        "-t", "--template_sdf", help="Template SDF.", default="template.sdf"
    )

    parser.add_argument("-o", "--output_sdf", help="Output SDF.", default="output.sdf")

    return parser.parse_args()


def fix_mol_conformer(rdkit_lig, rdkit_template, nthreads=1):

    mcsResult = rdFMCS.FindMCS(
        [rdkit_lig, rdkit_template],
        completeRingsOnly=True,
    )
    if mcsResult.smartsString and len(mcsResult.smartsString) > 0:
        patt = Chem.MolFromSmarts(mcsResult.smartsString, mergeHs=True)

        # keep only the MCS part of the reference molecule
        ref = AllChem.ReplaceSidechains(rdkit_template, patt)
        if ref:
            core = AllChem.DeleteSubstructs(ref, Chem.MolFromSmiles("*"))
            core.UpdatePropertyCache()

    rdkit_lig = Chem.AddHs(rdkit_lig)

    AllChem.ConstrainedEmbed(
        rdkit_lig, core
    )  # constrained minimization of newmol versus the core of the reference

    return rdkit_lig


def main():
    args = parse_args()

    rdkit_template = Chem.MolFromMolFile(args.template_sdf)
    if args.ligand_format == "smi":
        rdkit_lig_entries = Chem.rdmolfiles.SmilesMolSupplier(args.ligand)
    if args.ligand_format == "sdf":
        rdkit_lig_entries = Chem.rdmolfiles.SDMolSupplier(args.ligand)
    with Chem.SDWriter(args.output_sdf) as w:
        for rdkit_lig in rdkit_lig_entries:
            conf = fix_mol_conformer(rdkit_lig, rdkit_template)
            w.write(conf)

    return


if __name__ == "__main__":
    main()
