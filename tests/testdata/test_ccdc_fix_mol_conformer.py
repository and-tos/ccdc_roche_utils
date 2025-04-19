from pathlib import Path

import pytest
from ccdc import io

from ccdc_roche_utils import ccdc_fix_mol_conformer


@pytest.mark.parametrize(
    "ligand_sdf, template_sdf, expected_confs",
    [
        ("ligand_1.sdf", "template_1.sdf", 2),
        ("ligand_2.sdf", "template_2.sdf", 3),
        ("ligand_3.sdf", "template_3.sdf", 2),
    ],
)
def test_ccdc_fix_mol_conformer(ligand_sdf, template_sdf, expected_confs):

    template_sdf = str(Path(__file__).parent / template_sdf)
    ligand_sdf = str(Path(__file__).parent / ligand_sdf)

    ccdc_template = io.MoleculeReader(template_sdf)[0]
    ccdc_lig_entries = io.EntryReader(ligand_sdf)
    for ccdc_lig_entry in ccdc_lig_entries:
        conf_entries = ccdc_fix_mol_conformer.fix_mol_conformer(
            ccdc_lig_entry, ccdc_template, max_conformers=expected_confs
        )
        assert len(conf_entries) == expected_confs


def main():
    for args in [
        ("ligand_1.sdf", "template_1.sdf", 2),
        ("ligand_2.sdf", "template_2.sdf", 3),
        ("ligand_3.sdf", "template_3.sdf", 2),
    ]:
        test_ccdc_fix_mol_conformer(*args)


if __name__ == "__main__":
    main()
