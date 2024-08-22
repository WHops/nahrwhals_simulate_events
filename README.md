# Sequence Mutation and Generation

This repository contains tools to generate and mutate DNA sequences based on user-defined parameters. It includes scripts for simulation of sequences with segmental duplications, and for simulated chains of mutation based on NAHR between these repeats. The code here was developed to help benchmark the NAHRwhals tool. 

## Features

- **Sequence Generation:** Create sequences with segmental duplications of specific lengths and similarity percentages.
- **Sequence Mutation:** Introduce random NAHR-based mutations into sequences, supporting various mutation types like duplication, deletion, and inversion.

## Installation

Clone the repository and install the required packages:

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
pip install -r requirements.txt
```

## Usage
### Generate and Mutate Sequences

1. Generate Sequences:

Use generate_sim_sequences.py to generate sequences with specific lengths and similarity percentages:

```
python generate_sim_sequences.py <n_sequences> <output_folder> --segment_lengths <lengths> --segment_similarities <similarities>
```

Example:

```
python generate_sim_sequences.py 50 out --segment_lengths 100 500 --segment_similarities 0.9 0.95
```

2. Mutate Sequences:

Use mutate_sequences.py to introduce mutations into the generated sequences:

```
python mutate_sequences.py <sequence_file> <outputdir>
```

Example:

```
python mutate_sequences.py out/seq1_mode1_00_sim90_len500.fa mut
```

3. Automate the Process:

Run the generate_mutate.sh script to generate sequences and apply mutations automatically:

```
bash generate_mutate.sh
```


### Directory Structure
data/out/: Stores generated sequences.
data/mut/: Stores mutated sequences.