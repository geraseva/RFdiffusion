
# RFDiffusion — NA-Binding Protein Design

This is a fork of [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion) modified for the design of nucleic acid-binding proteins.

## Key modifications

1. **Nucleic acid support** — added handling for DNA/RNA chains with dedicated auxiliary potentials based on the `substrate_contacts` framework.

2. **Rigid-body averaging** — potentials computed per residue can now be averaged across the protein chain and applied as collective impulses for translation and rotation of the entire chain as a rigid body.

3. **Px0-based potential calculation** — instead of computing potentials on the noisy intermediate structure (`xt`), potentials are evaluated on `px0` (the extrapolated structure used for self-conditioning). Since `px0` is closer to a natural structure, this provides more accurate protein-NA interaction estimates.

## Usage

```bash
# Basic design with NA-binding auxiliary potentials
python ./scripts/run_inference.py \
    inference.output_prefix=examples/example_outputs/design_dsdna \
    inference.input_pdb=input_pdbs/5EDN.pdb \
    "contigmap.contigs=[70-130]" \
    potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth=0.6,predicted:1"] \
    potentials.guide_scale=1
```

For full documentation, see the [original RFDiffusion README](https://github.com/RosettaCommons/RFdiffusion).


## References

- Original RFDiffusion: [Watson et al., Nature 2023](https://doi.org/10.1038/s41586-023-06415-8)
- This fork: [Geraseva et al., JCIM 2026](https://doi.org/10.1021/acs.jcim.6c00109)
