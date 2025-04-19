#!/bin/bash

conda deactivate
conda activate ccdc_roche_3.9

ccdc_roche_3.9/bin/python /ccdc_roche_devel/modeling_utils/ccdc_fix_mol_conformer.py -l /ccdc_roche_devel/modeling_utils/tests/testdata/ligand_3.sdf -t /ccdc_roche_devel/modeling_utils/tests/testdata/template_3.sdf -o output.sdf
