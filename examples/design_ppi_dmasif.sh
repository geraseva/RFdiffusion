#!/bin/bash

../scripts/run_inference.py inference.output_prefix=example_outputs/design_ppi_dmasif_single \
  inference.input_pdb=input_pdbs/insulin_target.pdb 'contigmap.contigs=[70-100]' inference.num_designs=10 \
  'potentials.guiding_potentials=["type:dmasif_interactions,non_int_weight:0.5"]'