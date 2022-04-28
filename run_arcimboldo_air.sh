#!/bin/bash

python ./main.py \
--output_dir=/cri4/albert/Desktop/atzr_gio_fixed_pocket \
--fasta_path=/cri4/albert/Desktop/atzr_gio_fixed_pocket/atzr_tetramer.fasta \
--feature_mode=from_custom_pdb \
--custom_pdb_path=/cri4/albert/Desktop/atzr_gio_fixed_pocket/2_atzr_gio_RefinedTemplate.pdb \
--poly_ala_list_of_res_ranges=1-101,102-129,130-200,202-241,242-9000,

101-102,129-130,200-202,241-242,225-226,148-149
