import random
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import re

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def calculate_similarity(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / len(seq1)

def find_duplicated_regions(sequence, SD_len, similarity_threshold):
    seq_len = len(sequence)
    duplicated_regions = []

    goalposts = list(range(0, seq_len - SD_len, SD_len))
    for i in range(len(goalposts)):
        start1 = goalposts[i]
        window1 = sequence[start1:start1 + SD_len]

        for j in range(i + 1, len(goalposts)):
            start2 = goalposts[j]
            window2 = sequence[start2:start2 + SD_len]

            similarity_direct = calculate_similarity(window1, window2)
            if similarity_direct >= similarity_threshold:
                duplicated_regions.append((start1, start2, SD_len, similarity_direct, 'direct'))
                next
            
            reverse_window2 = reverse_complement(window2)
            similarity_reverse = calculate_similarity(window1, reverse_window2)
            if similarity_reverse >= similarity_threshold:
                duplicated_regions.append((start1, start2, SD_len, similarity_reverse, 'reverse'))
    
    return duplicated_regions

def mutate_sequence(sequence, sdsize, similarity_threshold, n_mutations, sequence_name, outputdir, seed=None):
    if seed is not None:
        random.seed(seed)

    sequence_list = [sequence]
    mutation_log = []
    
    for _ in range(n_mutations):
        duplicated_regions = find_duplicated_regions(sequence, sdsize, similarity_threshold)
        
        if not duplicated_regions:
            print("No more regions to mutate")
            return sequence, mutation_log
        
        for attempt in range(10):
            region = random.choice(duplicated_regions)
            start1, start2, length, similarity, orientation = region

            if orientation == 'direct':
                action = random.choice(["delete", "duplicate"])
            elif orientation == 'reverse':
                action = "invert"

            if action == "delete":
                sequence_candidate = delete_between(sequence, start1, start2)
            elif action == "duplicate":
                sequence_candidate = duplicate_between(sequence, start1, start2)
            elif action == 'invert':
                sequence_candidate = invert_between(sequence, start1, start2 + sdsize)

            if sequence_candidate in sequence_list:
                print("Sequence already in list. Retrying...")
                continue

            sequence = sequence_candidate
            mutation_log.append(f"{action}: {start1} to {start2}")
            sequence_list.append(sequence)

            sequence_name = sequence_name + f"_{action}_{start1}_{start2}"

            with open(f"{outputdir}/{sequence_name}.fa", "w") as output_handle:
                output_handle.write(f">{sequence_name}\n{sequence}")

            break

    return sequence, mutation_log

def delete_between(sequence, start, end):
    return sequence[:start] + sequence[end:]

def duplicate_between(sequence, start, end):
    to_duplicate = sequence[start:end]
    return sequence[:start] + to_duplicate + sequence[start:]

def invert_between(sequence, start, end):
    to_invert = sequence[start:end]
    inverted = reverse_complement(to_invert)
    return sequence[:start] + inverted + sequence[end:]

parser = argparse.ArgumentParser(description='Simulate mutations in a sequence.')
parser.add_argument('sequence_file', type=str, help='Path to the input FASTA file')
parser.add_argument('outputdir', type=str, help='Path to the output directory')
parser.add_argument('--seed', type=int, default=None, help='Seed for random number generator to ensure reproducibility.')

args = parser.parse_args()

sequence_file = args.sequence_file
outputdir = args.outputdir

pattern = r'.*/seq(\d+)_mode\d+_\d+_sim(\d+)_len(\d+)\.fa'

match = re.match(pattern, sequence_file)

if match:
    similarity_threshold = int(match.group(2)) / 100
    sdsize = int(int(match.group(3)))
else:
    raise ValueError("Filename does not match the expected format.")

for record in SeqIO.parse(sequence_file, "fasta"):
    sequence = str(record.seq)
    break

n_mutations = 3

sequence_name = sequence_file.split("/")[1].split(".")[0]

mutated_sequence, mutation_log = mutate_sequence(sequence, sdsize, similarity_threshold, n_mutations, sequence_name, outputdir, args.seed)

for log in mutation_log:
    print(log)