#!/bin/bash

# This is a benchmark file to reproduce experiments

source /home/domain/data/prog/miniconda3/etc/profile.d/conda.sh
RFDIFFUSION_PATH=../

conda activate rfdiffusion

# define rfdiffusion parameters

declare -A pot_dict
pot_dict[no_potential]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1"] potentials.guide_scale=0'
pot_dict[na_contacts]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1"] potentials.guide_scale=1'

pot_dict[na_contacts_02avg]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth:0.2"] potentials.guide_scale=1'
pot_dict[na_contacts_04avg]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth:0.4"] potentials.guide_scale=1'
pot_dict[na_contacts_06avg]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth:0.6"] potentials.guide_scale=1'
pot_dict[na_contacts_08avg]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth:0.8"] potentials.guide_scale=1'
pot_dict[na_contacts_avg]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,smooth:1"] potentials.guide_scale=1'

pot_dict[na_contacts_px0]='potentials.guiding_potentials=["type:na_contacts,s:2,r_0:8,rep_r_0:7,rep_s:2,d_0:2,rep_r_min:1,predicted:1"] potentials.guide_scale=1'

declare -A motif_dict
motif_dict[de_novo]="contigmap.contigs=[123-183]" 
motif_dict[site]="contigmap.contigs=[40-60/A21-22/40-60/B222-222/40-60] inference.ckpt_override_path=$RFDIFFUSION_PATH/models/ActiveSite_ckpt.pt" 
motif_dict[binding]="contigmap.contigs=[18-18/A24-24/18-18/A43-43/27-27/A71-71/3-3/A75-75/5-5/A81-81/85-85] inference.ckpt_override_path=$RFDIFFUSION_PATH/models/ActiveSite_ckpt.pt" 

mkdir example_outputs/na/
mkdir example_outputs/na/rfdiffusion_outputs

# run rfdiffusion
for motif in de_novo site binding;
do
  mkdir example_outputs/na/rfdiffusion_outputs/$motif
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    echo imsoi $motif $pot >> logs/test_na.log
    mkdir example_outputs/na/rfdiffusion_outputs/$motif/$pot
    $RFDIFFUSION_PATH/scripts/run_inference.py inference.num_designs=50 inference.deterministic=true \
    inference.output_prefix=example_outputs/na/rfdiffusion_outputs/$motif/$pot/imsoi \
    "inference.input_pdb=input_pdbs/1M5X.pdb" \
    ${motif_dict[$motif]} ${pot_dict[$pot]} >> example_outputs/na/rfdiffusion_outputs/test_na.log
  done 
done

conda deactivate

# add NA structures to designs

for motif in de_novo site binding;
do
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    for lig in `ls example_outputs/na/rfdiffusion_outputs/$motif/$pot| grep na`
    do
      cat example_outputs/na/rfdiffusion_outputs/$motif/$pot/$(echo $lig | sed 's/_na//') | egrep -v 'DA|DG|DC|DT' > tmp
      cat tmp > example_outputs/na/rfdiffusion_outputs/$motif/$pot/$(echo $lig | sed 's/_na//')
    numatoms=$(tail -1 example_outputs/na/rfdiffusion_outputs/$motif/$pot/$(echo $lig | sed 's/_na//') | head -c 11 | tail -c 4 | xargs)
    cat example_outputs/na/rfdiffusion_outputs/$motif/$pot/$lig | \
      awk -v var=$numatoms 'BEGIN { FIELDWIDTHS="6 5 2 1 54"}{printf "%s%5u%s%s%s%s          %s\n", $1,$2 + var,$3,$4,$5,$6,$4}' \
      >> example_outputs/na/rfdiffusion_outputs/$motif/$pot/$(echo $lig | sed 's/_na//')
    done
  done
done 


# running LigandMPNN
LIGANDMPNN_PATH=../../LigandMPNN/
WD=`pwd`

conda activate ligandmpnn_env

mkdir example_outputs/na/proteinmpnn_output

for motif in de_novo site binding;
do
  mkdir example_outputs/na/proteinmpnn_output/$motif

  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    # generate json files for ligandmpnn multiple run
    cd $WD
    pdb_ids=$WD/example_outputs/na/proteinmpnn_output/$motif/${pot}_ids.json
    redesigned_residues=$WD/example_outputs/na/proteinmpnn_output/$motif/${pot}_res.json
    echo { > ${pdb_ids}
    echo { > ${redesigned_residues}

    itdir=$WD/example_outputs/na/proteinmpnn_output/$motif/$pot
    for pdb in `ls $itdir | grep pdb | grep -v na `;
    do
      echo \"$itdir/$pdb\": \"\", >> ${pdb_ids}
      echo \"$itdir/$pdb\": \"$(cat $itdir/$pdb | grep 'CA  GLY A' | grep '0.00$' | awk 'BEGIN{OFS="";ORS=" "}{print $5, $6}' | sed 's/\ $//' )\", >> ${redesigned_residues}
    done
  
    tmp=`cat ${pdb_ids} | sed '$ s/\,$/\n}/'`
    echo $tmp > ${pdb_ids}
    tmp=`cat ${redesigned_residues} | sed '$ s/\,$/\n}/'`
    echo $tmp > ${redesigned_residues}

    # run ligandmpnn
    cd $LIGANDMPNN_PATH
    python3 run.py \
    --pdb_path_multi ${pdb_ids} \
    --out_folder $WD/example_outputs/na/proteinmpnn_output/$motif/$pot \
    --model_type "ligand_mpnn" \
    --batch_size 5 --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1
    --redesigned_residues_multi ${redesigned_residues}

   # write ligandmpnn output sequences to a single file
    cd $WD
    for pdb in `ls example_outputs/na/proteinmpnn_output/$motif/$pot/seqs | grep fa` ;
    do
      echo $pdb
      cat example_outputs/na/proteinmpnn_output/$motif/$pot/seqs/${pdb} | awk '{print $1, $2}' | grep id= -A 1 | sed 's/,\ id=/_/' | sed 's/,//' >> example_outputs/na/proteinmpnn_output/$motif/$pot/all_seqs.fasta   
    done
  done
done
conda deactivate

# relax designs
conda activate pyrosetta

for motif in de_novo site binding;
do
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    mkdir example_outputs/na/proteinmpnn_output/${motif}/${pot}/relaxed
    python $RFDIFFUSION_PATH/scripts/fast_relax.py -pdbdir example_outputs/na/proteinmpnn_output/${motif}/${pot}/backbones \
    -outdir example_outputs/na/proteinmpnn_output/${motif}/${pot}/relaxed -parallel 
  done
done

conda deactivate

# calculate secondary structures

DSSP_BIN=/home/domain/anur/progs/bin/dsspcmbi
mkdir example_outputs/na/dssp_statistics

for motif in de_novo site binding;
do
  mkdir example_outputs/na/dssp_statistics/$motif/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    mkdir example_outputs/na/dssp_statistics/$motif//$pot
    for f in `ls example_outputs/na/proteinmpnn_output/$motif/$pot/relaxed/| grep pdb`
    do 
      $DSSP_BIN example_outputs/na/proteinmpnn_output/$motif/$pot/relaxed/$f example_outputs/na/dssp_statistics/$motif//$pot/$f.dssp
    done
  done
done


# calculate NA binding surfaces

conda activate dmasif
DMASIF_PATH=../../masif_npi

mkdir example_outputs/na/dmasif_output/

for motif in de_novo site binding;
do
  mkdir example_outputs/na/dmasif_output/$motif/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0; 
  do
    mkdir example_outputs/na/dmasif_output/$motif/$pot
    ls example_outputs/na/proteinmpnn_output/$motif/$pot/relaxed/ | awk '{print $1,"A","CD"}' > example_outputs/na/dmasif_output/$motif/$pot.list
    python3 $DMASIF_PATH/train_inf.py inference --experiment_name npi_site_b2 --site  --na NA --protonate \
      --pdb_list example_outputs/na/dmasif_output/$motif/$pot.list --data_dir example_outputs/na/proteinmpnn_output/$motif/$pot/relaxed/ \
      --out_dir example_outputs/na/dmasif_output/$motif/$pot
    done
  done
done 

END_COMMENT


# fold complexes using boltz2

conda activate boltz2

mkdir example_outputs/na/boltz_predictions

for motif in de_novo;
do
  mkdir example_outputs/na/boltz_predictions/$motif/
  for pot in na_contacts_px0; 
  do
    mkdir example_outputs/na/boltz_predictions/$motif/$pot
    mkdir example_outputs/na/boltz_predictions/$pot/conf
    
    cat example_outputs/na/proteinmpnn_output/$motif/$pot/all_seqs.fasta | while read name; read seq
    do 
      echo ${name:1}
      cat > example_outputs/na/boltz_predictions/$motif/$pot/conf/${name:1}.yaml << EOL
sequences:
  - dna:
      id: C
      sequence: GCAGAACGTCGTGAGACAGTTCCG
  - dna:
      id: D
      sequence: CGGAACTGTCTCACGACGTTCTGC
  - protein:
      id: A
      msa: empty
EOL
    echo '      sequence:' ${seq} >> example_outputs/na/boltz_predictions/$motif/$pot/conf/${name:1}.yaml 
    done
  boltz predict example_outputs/na/boltz_predictions/$motif/$pot/conf/ --out_dir example_outputs/na/boltz_predictions/$motif/$pot/ --output_format pdb
  done
done

conda deactivate
