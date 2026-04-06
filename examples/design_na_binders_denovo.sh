#!/bin/bash

# This is a benchmark file to reproduce experiments

source /home/domain/data/prog/miniconda3/etc/profile.d/conda.sh
eval "$(/home/domain/data/prog/micromamba/bin/micromamba shell hook -s posix)"

RFDIFFUSION_PATH=../
WD=`pwd`

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

declare -A pdb_dict
pdb_dict[dsDNA-IMsoI]=input_pdbs/1M5X.pdb
pdb_dict[dsDNA-HoxB13]=input_pdbs/5EDN.pdb 
pdb_dict[DNAapt-LDH]=input_pdbs/5HTO.pdb 
pdb_dict[gqDNA-NB]=input_pdbs/9GXH.pdb 
pdb_dict[TATA-TBP]=input_pdbs/1QNA.pdb 
pdb_dict[ssDNA-Enc34]=input_pdbs/5ODL.pdb 

pdb_dict[RNAapt-Thrb]=input_pdbs/5DO4.pdb 
pdb_dict[miRNA-LIN28A]=input_pdbs/5UDZ.pdb 
pdb_dict[dsRNA-hStau1]=input_pdbs/6HTU.pdb 
pdb_dict[snRNA-U2]=input_pdbs/1A9N.pdb 
pdb_dict[RNAapt-IgG]=input_pdbs/3AGV.pdb 

declare -A motif_dict
motif_dict[dsDNA-IMsoI]="contigmap.contigs=[6-36/A24-24/18-18/A43-43/27-27/A71-71/3-3/A75-75/5-5/A81-81/6-36]"
motif_dict[dsDNA-HoxB13]="contigmap.contigs=[29-59/A258-258/3-3/A262-262/6-6/A269-269/29-59]"
motif_dict[DNAapt-LDH]="contigmap.contigs=[3-33/B14-14/25-25/B40-40/3-3/B44-44/30-30/B82-84/3-33]"
motif_dict[gqDNA-NB]="contigmap.contigs=[6-36/B37-37/1-1/B39-39/4-4/B44-45/1-1/B47-47/35-35/B94-94/9-9/B104-105/6-36]"
motif_dict[TATA-TBP]="contigmap.contigs=[7-37/A27-27/1-1/A29-29/27-27/A57-57/14-14/A72-72/1-1/A74-74/5-5/A80-80/1-1/A82-82/7-37]"
motif_dict[ssDNA-Enc34]="contigmap.contigs=[3-33/A43-43/30-30/A112-112/27-27/A150-150/3-3/A154-154/3-33]"

motif_dict[RNAapt-Thrb]="contigmap.contigs=[4-34/H98-98/25-25/H169-169/4-4/H174-175/25-25/H245-245/2-2/H248-248/4-34]"
motif_dict[miRNA-LIN28A]="contigmap.contigs=[5-35/A45-46/8-8/A55-55/15-15/A71-71/3-3/A75-75/2-2/A78-78/5-5/A84-85/14-14/A100-100/1-1/A102-102/2-2/A105-105/4-34]"
motif_dict[dsRNA-hStau1]="contigmap.contigs=[6-36/A183-184/2-2/A187-187/3-3/A191-191/3-3/A195-195/15-15/A211-212/21-21/A234-235/2-2/A238-239/7-37]"
motif_dict[snRNA-U2]="contigmap.contigs=[19-49/B19-19/24-24/B44-49/20-50]"
motif_dict[RNAapt-IgG]="contigmap.contigs=[4-34/B340-340/1-1/B342-342/1-1/B344-344/28-28/B373-374/27-27/B402-402/3-33]"



conda activate rfdiffusion

mkdir example_outputs/na/
mkdir example_outputs/na/rfdiffusion_outputs
mkdir example_outputs/na/logs

# run rfdiffusion
for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg na_contacts_px0 na_contacts_px0_06avg; 
  do
    echo "Start ${complex} $pot" >> example_outputs/na/logs/rfdiffusion.log
    mkdir example_outputs/na/rfdiffusion_outputs/${complex}/$pot
    $RFDIFFUSION_PATH/scripts/run_inference.py inference.num_designs=50 inference.deterministic=true \
    inference.output_prefix=example_outputs/na/rfdiffusion_outputs/${complex}/$pot/design \
    inference.input_pdb=${pdb_dict[${complex}]} \
    "contigmap.contigs=[70-130]" ${pot_dict[$pot]} >> example_outputs/na/logs/rfdiffusion.log

    # make complexes
    mkdir example_outputs/na/rfdiffusion_outputs/${complex}/$pot/apo
    mkdir example_outputs/na/rfdiffusion_outputs/${complex}/$pot/complex
    mkdir example_outputs/na/rfdiffusion_outputs/${complex}/$pot/na
    for lig in `ls example_outputs/na/rfdiffusion_outputs/${complex}/$pot| grep _na`
    do
    cp example_outputs/na/rfdiffusion_outputs/${complex}/$pot/$(echo $lig | sed 's/_na//') \
      example_outputs/na/rfdiffusion_outputs/${complex}/$pot/complex/$(echo $lig | sed 's/_na//')
    numatoms=$(tail -1 example_outputs/na/rfdiffusion_outputs/${complex}/$pot/$(echo $lig | sed 's/_na//') | head -c 11 | tail -c 4 | xargs)
    cat example_outputs/na/rfdiffusion_outputs/${complex}/$pot/$lig | \
      awk -v var=$numatoms 'BEGIN { FIELDWIDTHS="6 5 2 1 54"}{printf "%s%5u%s%s%s%s          %s\n", $1,$2 + var,$3,$4,$5,$6,$4}' \
      >> example_outputs/na/rfdiffusion_outputs/${complex}/$pot/complex/$(echo $lig | sed 's/_na//')
    mv example_outputs/na/rfdiffusion_outputs/${complex}/$pot/$(echo $lig | sed 's/_na//') \
      example_outputs/na/rfdiffusion_outputs/${complex}/$pot/apo/$(echo $lig | sed 's/_na//')
    mv example_outputs/na/rfdiffusion_outputs/${complex}/$pot/$lig \
      example_outputs/na/rfdiffusion_outputs/${complex}/$pot/na/$(echo $lig | sed 's/_na//')    
    done
  done 
done

# run motif scaffolding

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}
  for pot in no_potential na_contacts; 
  do
    echo "Start ${complex} mk_$pot" >> example_outputs/na/logs/rfdiffusion_mk.log
    mkdir example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot
    $RFDIFFUSION_PATH/scripts/run_inference.py inference.num_designs=50 inference.deterministic=true \
    inference.output_prefix=example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot/design \
    inference.input_pdb=${pdb_dict[${complex}]} ${pot_dict[$pot]} \
    ${motif_dict[$complex]} inference.ckpt_override_path=$RFDIFFUSION_PATH/models/ActiveSite_ckpt.pt >> example_outputs/na/logs/rfdiffusion_mk.log

    # add NA structures to designs
    for lig in `ls example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot| grep _na`
    do
    numatoms=$(tail -1 example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot/$(echo $lig | sed 's/_na//') | head -c 11 | tail -c 4 | xargs)
    cat example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot/$lig | \
      awk -v var=$numatoms 'BEGIN { FIELDWIDTHS="6 5 2 1 54"}{printf "%s%5u%s%s%s%s          %s\n", $1,$2 + var,$3,$4,$5,$6,$4}' \
      >> example_outputs/na/rfdiffusion_outputs/${complex}/mk_$pot/$(echo $lig | sed 's/_na//')
    done
  done 
done

conda deactivate


# run rfdiffusion3

micromamba activate /mnt/storage/prog/micromamba/envs/foundry

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}/rfd3
   
  contig=$(cat ${pdb_dict[${complex}]} | grep "C1'" | awk '{
    chain = substr($0,22,1)
    resnum = substr($0,23,4) + 0
    
    if (!min[chain] || resnum < min[chain]) min[chain] = resnum
    if (!max[chain] || resnum > max[chain]) max[chain] = resnum
    } END {
    contig = ""
    for (c in min) {
        if (contig != "") contig = contig ",/0,"
        contig = contig c min[c] "-" max[c]
    }
    print contig
    }')
  length=$(cat ${pdb_dict[${complex}]} | grep "C1'" | wc -l)
  echo '{"design": {' > example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  echo '  "input": "'$WD/${pdb_dict[${complex}]}'",' >> example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  echo '  "contig": "70-130,/0,'${contig}'",' >> example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  echo '  "length": "'$(($length + 70))-$(($length + 130))'"' >> example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  echo '}}' >> example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  cat example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json

  rfd3 design out_dir=example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/ \
     diffusion_batch_size=1 n_batches=50 skip_existing=True \
     inputs=example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/config.json
  
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/apo
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/complex
  mkdir example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/na

  for file in `ls example_outputs/na/rfdiffusion_outputs/${complex}/rfd3 | grep cif.gz `;
    do
      gzip -d example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/$file
    done
  for file in `ls example_outputs/na/rfdiffusion_outputs/${complex}/rfd3 | grep cif `;
    do
      python /home/domain/geraseva/geraseva_data/auxiliary_potential/utils/cif2pdb.py \
         example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/$file \
         $(echo example_outputs/na/rfdiffusion_outputs/${complex}/rfd3/$file | sed 's/...$/pdb/')
      cat example_outputs/na/rfdiffusion_outputs/$complex/rfd3/$(echo $file | sed 's/...$/pdb/') | sed -n '/^TER/q; p' | \
        grep -E 'ATOM|HETATM' > example_outputs/na/rfdiffusion_outputs/$complex/rfd3/apo/$(echo $file | sed 's/...$/pdb/')
      cat example_outputs/na/rfdiffusion_outputs/$complex/rfd3/$(echo $file | sed 's/...$/pdb/') | sed -n '/^TER/,$p' | \
        grep -E 'ATOM|HETATM' > example_outputs/na/rfdiffusion_outputs/$complex/rfd3/na/$(echo $file | sed 's/...$/pdb/')

    done
done
micromamba deactivate

# running LigandMPNN

LIGANDMPNN_PATH=../../LigandMPNN/
WD=`pwd`

conda activate ligandmpnn_env

mkdir example_outputs/na/proteinmpnn_output

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/proteinmpnn_output/${complex}

  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
  do

    # generate json files for ligandmpnn multiple run
    cd $WD
    pdb_ids=$WD/example_outputs/na/proteinmpnn_output/${complex}/${pot}_ids.json
    redesigned_residues=$WD/example_outputs/na/proteinmpnn_output/${complex}/${pot}_res.json
    echo { > ${pdb_ids}
    echo { > ${redesigned_residues}

    itdir=$WD/example_outputs/na/rfdiffusion_outputs/${complex}/$pot
    for pdb in `ls $itdir | grep pdb | grep -v na `;
    do
      echo \"$itdir/$pdb\": \"\", >> ${pdb_ids}
      echo \"$itdir/$pdb\": \"$(cat $itdir/$pdb | grep '^ATOM.\{9\}CA' | grep -v '1.00$' | awk 'BEGIN{OFS="";ORS=" "}{print $5, $6}' | sed 's/\ $//' )\", >> ${redesigned_residues}
    done
  
    tmp=`cat ${pdb_ids} | sed '$ s/\,$/\n}/'` 
    echo $tmp > ${pdb_ids}
    tmp=`cat ${redesigned_residues} | sed '$ s/\,$/\n}/'`
    echo $tmp  > ${redesigned_residues}

    # run ligandmpnn
    cd $LIGANDMPNN_PATH
    python3 run.py \
    --pdb_path_multi ${pdb_ids} \
    --out_folder $WD/example_outputs/na/proteinmpnn_output/${complex}/${pot} \
    --model_type "ligand_mpnn" \
    --batch_size 5 --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
    --redesigned_residues_multi ${redesigned_residues}

   # write ligandmpnn output sequences to a single file
    cd $WD
    echo -n > example_outputs/na/proteinmpnn_output/${complex}/${pot}/all_seqs.fasta   
    for pdb in `ls example_outputs/na/proteinmpnn_output/${complex}/${pot}/seqs | grep fa` ;
    do
      echo $pdb
      cat example_outputs/na/proteinmpnn_output/${complex}/${pot}/seqs/${pdb} | awk '{print $1, $2}' | grep id= -A 1 | sed 's/,\ id=/_/' | sed 's/,//' >> example_outputs/na/proteinmpnn_output/${complex}/$pot/all_seqs.fasta   
    done
  done
done
conda deactivate


# relax designs
conda activate pyrosetta

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
  do
    mkdir example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed
    python $RFDIFFUSION_PATH/scripts/fast_relax.py -pdbdir example_outputs/na/proteinmpnn_output/${complex}/${pot}/backbones \
    -outdir example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed -parallel 
  done
done

conda deactivate

# calculate secondary structures

DSSP_BIN=/home/domain/anur/progs/bin/dssp
mkdir example_outputs/na/dssp_statistics

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/dssp_statistics/$complex/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
    mkdir example_outputs/na/dssp_statistics/$complex/$pot
    for f in `ls example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/| grep pdb`
    do 
      $DSSP_BIN example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/$f example_outputs/na/dssp_statistics/$complex/$pot/$f.dssp
    done
  done
done

# calculate structural diversity

micromamba activate /mnt/storage/prog/micromamba/envs/foldseek

mkdir example_outputs/na/foldseek_clustering
mkdir example_outputs/na/foldseek_clustering/tmp
#foldseek databases PDB /home/domain/data/geraseva/pdb/pdb \
#  example_outputs/na/foldseek_clustering/tmp

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/foldseek_clustering/$complex/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
  do
    mkdir example_outputs/na/foldseek_clustering/$complex/$pot
    foldseek easy-search example_outputs/na/rfdiffusion_outputs/$complex/$pot/apo \
      /home/domain/data/geraseva/pdb/pdb --alignment-type 0 \
      --format-output "query,target,alntmscore,qtmscore,lddt,prob" \
      example_outputs/na/foldseek_clustering/$complex/$pot/design \
      example_outputs/na/foldseek_clustering/tmp/ 
  done
done

micromamba deactivate

# calculate ddg


for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg na_contacts_avg \
    na_contacts_px0 na_contacts_px0_06avg mk_no_potential mk_na_contacts rfd3;
  do
    mkdir example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed_apo/
    mkdir example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed_na/
    for file in `ls example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed/ | grep pdb`;
    do
      cat example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed/${file} | \
        sed -n '/^TER/q; p' | grep -E 'ATOM|HETATM' > example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed_apo/${file}
      cat example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed/${file} | \
        sed -n '/^TER/,$p' | grep -E 'ATOM|HETATM' | sed 's/./C/22' > example_outputs/na/proteinmpnn_output/${complex}/${pot}/relaxed_na/${file}
    done
  done
done
COMMENT

conda activate pyrosetta

python3 - << 'END_SCRIPT'
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

from pyrosetta import *
from pyrosetta.rosetta import *
#from pyrosetta.toolbox import *
from pyrosetta.teaching import *

init("-mute all")

scorefxnDDG=get_fa_scorefxn()
fa_sol_score=ScoreFunction()
fa_sol_score.set_weight(fa_sol, 1.0)
interface_score=ScoreFunction()
interface_score.set_weight(interface_dd_pair, 1.0)

from multiprocessing import Pool, cpu_count
from IPython.display import clear_output


def calc_ddg(fs):
    testPose= Pose()
    testPose = pose_from_pdb(fs)
    ddg=scorefxnDDG(testPose)
    clear_output()
    return ddg

na_pots=['no_potential','na_contacts','na_contacts_02avg','na_contacts_04avg','na_contacts_06avg','na_contacts_08avg','na_contacts_avg','na_contacts_px0',
         'na_contacts_px0_06avg','mk_no_potential','mk_na_contacts','rfd3']
complexes=['dsDNA-IMsoI', 'dsDNA-HoxB13', 'DNAapt-LDH', 'gqDNA-NB', 'TATA-TBP', 'ssDNA-Enc34', 'RNAapt-Thrb', 'miRNA-LIN28A', 'dsRNA-hStau1', 'snRNA-U2', 'RNAapt-IgG']

arr=[]
for complex in complexes:
    if not os.path.isdir(f'example_outputs/na/rfdiffusion_outputs/{complex}'):
        continue
    for potential in na_pots:
        if not os.path.isdir(f'example_outputs/na/proteinmpnn_output/{complex}/{potential}/relaxed/'):
            continue
        for file in os.listdir(f'example_outputs/na/proteinmpnn_output/{complex}/{potential}/relaxed/'):
            if file[-3:]!='pdb':
                continue
            arr.append([f'example_outputs/na/proteinmpnn_output/{complex}/{potential}/relaxed/{file}', 
                        f'example_outputs/na/proteinmpnn_output/{complex}/{potential}/relaxed_apo/{file}',
                        f'example_outputs/na/proteinmpnn_output/{complex}/{potential}/relaxed_na/{file}', 
                        file.split('.')[0],complex,potential])
        
        
str_df=pd.DataFrame(columns=['path_c','path_a','path_n', 'name','complex','potential'], 
                    data=arr)

with Pool(cpu_count()) as p:
    scores_c=np.array(list(tqdm(p.imap(calc_ddg, str_df.path_c), total=len(str_df.path_c))))
    scores_a=np.array(list(tqdm(p.imap(calc_ddg, str_df.path_a), total=len(str_df.path_a))))
    scores_n=np.array(list(tqdm(p.imap(calc_ddg, str_df.path_n), total=len(str_df.path_n))))

str_df['ddg']=scores_c-(scores_a+scores_n)
str_df.drop('path_c', axis=1, inplace=True)
str_df.drop('path_a', axis=1, inplace=True)
str_df.drop('path_n', axis=1, inplace=True)

str_df.to_csv('example_outputs/na/logs/str_df.csv',index=False)
END_SCRIPT

conda deactivate

# fold complexes using boltz2

conda activate boltz2

mkdir example_outputs/na/boltz_predictions

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/boltz_predictions/$complex/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
  do
    mkdir example_outputs/na/boltz_predictions/$complex/$pot
    mkdir example_outputs/na/boltz_predictions/$complex/$pot/conf
    
    cat example_outputs/na/proteinmpnn_output/$complex/$pot/all_seqs.fasta | while read name; read seq
    do 
      echo ${name:1}
      
      cat > example_outputs/na/boltz_predictions/$complex/$pot/conf/${name:1}.yaml << EOL
sequences:
  - protein:
      id: A
      msa: empty
EOL
    echo '      sequence:' ${seq} >> example_outputs/na/boltz_predictions/$complex/$pot/conf/${name:1}.yaml 
    cat ${pdb_dict[${complex}]} | grep "C1'" | \
        awk '{
        c=substr($0,22,1)
        b=substr($0,18,4)
        t=(b~/ D|^D/)?"dna":"rna"
        gsub(/[ D]/,"",b)
        s[c,t]=s[c,t] substr(b,1,1)
        } END {
        for(ct in s) {
          split(ct, a, SUBSEP)
          print "  - " a[2] ":\n      id: " a[1] "\n      sequence: " s[ct]
          }
        }' >> example_outputs/na/boltz_predictions/$complex/$pot/conf/${name:1}.yaml
    done
  boltz predict example_outputs/na/boltz_predictions/$complex/$pot/conf/ --out_dir example_outputs/na/boltz_predictions/$complex/$pot/ --output_format pdb
  done
done

conda deactivate

# calculate NA binding surfaces

conda activate dmasif
DMASIF_PATH=/home/domain/data/geraseva/masif_npi

mkdir example_outputs/na/dmasif_output/

for complex in dsDNA-IMsoI dsDNA-HoxB13 DNAapt-LDH gqDNA-NB TATA-TBP ssDNA-Enc34 RNAapt-Thrb miRNA-LIN28A dsRNA-hStau1 snRNA-U2 RNAapt-IgG;
do
  mkdir example_outputs/na/dmasif_output/$complex/
  for pot in no_potential na_contacts na_contacts_02avg na_contacts_04avg na_contacts_06avg na_contacts_08avg \
    na_contacts_avg na_contacts_px0 na_contacts_px0_06avg rfd3 mk_no_potential mk_na_contacts; 
  do
    mkdir example_outputs/na/dmasif_output/$complex/$pot
    protchains=$(cat example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/*_0_1.pdb | \
            grep '^ATOM.\{9\}CA' | grep -o '^ATOM.\{17\}\(.\)' | cut -c22 | sort -u | tr -d '\n')
    chains=$(cat example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/*_0_1.pdb | \
            grep "^ATOM.\{9\}C1'" | grep -o '^ATOM.\{17\}\(.\)' | cut -c22 | sort -u | tr -d '\n')           
    ls example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/ | \
      awk -v chains=$chains -v prot=$protchains '{print $1,prot,chains}' > example_outputs/na/dmasif_output/$complex/$pot.list
    python3 $DMASIF_PATH/train_inf.py inference --experiment_name npi_site_b2 --site  --na NA \
      --pdb_list example_outputs/na/dmasif_output/$complex/$pot.list \
      --data_dir example_outputs/na/proteinmpnn_output/$complex/$pot/relaxed/ \
      --out_dir example_outputs/na/dmasif_output/$complex/$pot   
  done
done 
conda deactivate



