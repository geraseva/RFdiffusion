#!/bin/bash

echo ${PYTHONPATH}

#../scripts/run_inference.py inference.output_prefix=example_outputs/design_ppi_dmasif_single \
#  inference.input_pdb=input_pdbs/insulin_target.pdb 'contigmap.contigs=[70-100]' inference.num_designs=10 \
#  'potentials.guiding_potentials=["type:dmasif_interactions,non_int_weight:0.5"]'

../scripts/run_inference.py inference.output_prefix=example_outputs/design_ppi_dmasif_dimer \
  inference.input_pdb=input_pdbs/insulin_target.pdb 'contigmap.contigs=[A1-150/0 70-100]' inference.num_designs=10 \
 'potentials.guiding_potentials=["type:dmasif_interactions,non_int_weight:50,int_weight:50"]' \
  potentials.guide_scale=2 potentials.guide_decay='quadratic'