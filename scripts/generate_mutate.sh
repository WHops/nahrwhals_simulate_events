#!/bin/bash

python generate_sim_sequences.py 50 out

mkdir -p mut
for f in out/*.fa
do 
    echo "Processing $f file.."; 
    python mutate_sequences.py "$f" mut
    python mutate_sequences.py "$f" mut
    python mutate_sequences.py "$f" mut
done