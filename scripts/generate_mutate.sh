#!/bin/bash

SEED=42  # Set a seed value for reproducibility

python generate_sim_sequences.py 1 out --seed $SEED

mkdir -p mut
for f in out/*.fa
do 
    echo "Processing $f file.."; 
    python mutate_sequences.py "$f" mut --seed $SEED
    python mutate_sequences.py "$f" mut --seed $SEED
    python mutate_sequences.py "$f" mut --seed $SEED
done
