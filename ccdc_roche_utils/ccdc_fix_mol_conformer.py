from ccdc import io, conformer, descriptors, molecule, entry
import argparse
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from ccdc_roche_scoring import template_docking_mcs
import time


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Generate conformers with fixed fragment.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument("-l", "--ligand_sdf", help="Ligand SDF.", default="ligand.sdf")

    parser.add_argument(
        "-t", "--template_sdf", help="Template SDF.", default="template.sdf"
    )

    parser.add_argument("-o", "--output_sdf", help="Output SDF.", default="output.sdf")

    return parser.parse_args()


def _rdkit_overlay(ccdc_lig_mol, ccdc_template):
    rdkit_lig = Chem.MolFromMol2Block(ccdc_lig_mol.to_string())
    rdkit_template = Chem.MolFromMol2Block(ccdc_template.to_string())

    mcsResult = rdFMCS.FindMCS(
        [rdkit_lig, rdkit_template], completeRingsOnly=True
    )  # find the maximum common substructure

    if mcsResult.smartsString and len(mcsResult.smartsString) > 0:
        patt = Chem.MolFromSmarts(mcsResult.smartsString, mergeHs=True)

        # keep only the core of the reference molecule
        ref = AllChem.ReplaceSidechains(rdkit_template, patt)
        if ref:
            core = AllChem.DeleteSubstructs(ref, Chem.MolFromSmiles("*"))
            core.UpdatePropertyCache()

    newmol = Chem.Mol(
        rdkit_lig
    )  # create a new instance of the ligand, as we will change the coordinates
    newmol = Chem.AddHs(newmol)
    AllChem.ConstrainedEmbed(
        newmol, core
    )  # constrained minimization of newmol versus the core of the reference

    ccdc_lig_mol = molecule.Molecule.from_string(
        Chem.MolToMolBlock(newmol), format="sdf"
    )
    return ccdc_lig_mol


def fix_mol_conformer(ccdc_lig_entry, ccdc_template, max_conformers=5, nthreads=1):

    ccdc_lig_mol = ccdc_lig_entry.molecule

    t1 = time.time()
    ccdc_lig_mol = _rdkit_overlay(ccdc_lig_mol, ccdc_template)
    t2 = time.time()
    print("RDKit overlay: ", t2 - t1)

    mcs = template_docking_mcs.MCS(ccdc_lig_mol, ccdc_template)
    template, mcs_atoms, mcs_bonds = mcs.return_mcs_scaffold(
        partial_ring_matches_allowed=False, ignore_hydrogens=True
    )
    overlayer = descriptors.MolecularDescriptors.Overlay(
        ccdc_template, ccdc_lig_mol, mcs_atoms
    )

    # lock MCS bonds
    ccdc_lig_mol = overlayer.molecule
    mcs = template_docking_mcs.MCS(ccdc_lig_mol, ccdc_template)
    template, mcs_atoms, mcs_bonds = mcs.return_mcs_scaffold(
        partial_ring_matches_allowed=False, ignore_hydrogens=True
    )
    lig_mcs_bonds = [bonds[1] for bonds in mcs_bonds]
    t3 = time.time()
    print("MCS matching: ", t3 - t2)

    conf_gen = conformer.ConformerGenerator(nthreads=nthreads)
    conf_gen.settings.max_conformers = max_conformers

    for lig_bond in lig_mcs_bonds:
        conf_gen.lock_torsion(lig_bond)

    confs = conf_gen.generate(ccdc_lig_mol)

    conf_entries = []
    for conf_cnt, conf in enumerate(confs):
        conf = conf.molecule
        conf.kekulize()
        conf = entry.Entry.from_molecule(conf)
        conf.attributes = ccdc_lig_entry.attributes
        conf.attributes["Conformer count"] = conf_cnt
        conf_entries.append(conf)
    t4 = time.time()
    print("Conf gen: ", t4 - t3)
    return conf_entries


def main():
    args = parse_args()
    ccdc_template = io.MoleculeReader(args.template_sdf)[0]
    ccdc_lig_entries = io.EntryReader(args.ligand_sdf)
    with io.EntryWriter(args.output_sdf) as w:
        for ccdc_lig_entry in ccdc_lig_entries:
            conf_entries = fix_mol_conformer(ccdc_lig_entry, ccdc_template)
            for conf in conf_entries:
                w.write(conf)

    return


if __name__ == "__main__":
    main()
