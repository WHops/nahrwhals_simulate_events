import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import nt_search
import pdb
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

    # Generate goalposts
    goalposts = list(range(0, seq_len - SD_len, SD_len))
    for i in range(len(goalposts)):
        start1 = goalposts[i]
        window1 = sequence[start1:start1 + SD_len]

        for j in range(i + 1, len(goalposts)):
            start2 = goalposts[j]
            window2 = sequence[start2:start2 + SD_len]
            
            # Check for direct match
            similarity_direct = calculate_similarity(window1, window2)
            if similarity_direct >= similarity_threshold:
                duplicated_regions.append((start1, start2, SD_len, similarity_direct, 'direct'))
                next
            
            # Check for reverse complement match
            reverse_window2 = reverse_complement(window2)
            similarity_reverse = calculate_similarity(window1, reverse_window2)
            if similarity_reverse >= similarity_threshold:
                duplicated_regions.append((start1, start2, SD_len, similarity_reverse, 'reverse'))
    
    return duplicated_regions


def mutate_sequence(sequence, sdsize, similarity_threshold, n_mutations, sequence_name, outputdir):

    # Make a list in which we store the sequence  
    sequence_list = [sequence]

    mutation_log = []
    
    for _ in range(n_mutations):

        duplicated_regions = find_duplicated_regions(sequence, sdsize, similarity_threshold)
        
        print(duplicated_regions)
        if not duplicated_regions:
            print("No more regions to mutate")
            return sequence, mutation_log  # No more regions to mutate
        
        for attempt in range(10):  # Retry up to 10 times if a conflicting mutation is detected

            region = random.choice(duplicated_regions)
            start1, start2, length, similarity, orientation = region

            if orientation=='direct':  # Same orientation
                action = random.choice(["delete", "duplicate"])
            elif orientation=='reverse':
                action = "invert"

            if action == "delete":
                sequence_candidate = delete_between(sequence, start1, start2)
            elif action == "duplicate":
                sequence_candidate = duplicate_between(sequence, start1, start2)
            elif action=='invert':  # Reverse orientation
                sequence_candidate = invert_between(sequence, start1, start2+sdsize)

            # test if the sequence is in sequence_list
            if sequence_candidate in sequence_list:
                print("Sequence already in list. Retrying...")
                continue

            # Infom the usr about the mutation
            print(f"Mutating region {start1} to {start2} with action {action}")


            sequence = sequence_candidate

            mutation_log.append(f"{action}: {start1} to {start2}")

            # Save the sequence to the list
            sequence_list.append(sequence)

            sequence_name = sequence_name + f"_{action}_{start1}_{start2}"
            # save the sequence as fasta using bioython, and using the sequence_name
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


# Sample sequence of length 20 with window size 5 for testing
#sequence = "AGCTAGCTAGAGCTTGATCGAGCTTGTAGCCGATCAAGCTTGATGCTGGA"

# Set up argument parser
parser = argparse.ArgumentParser(description='Simulate mutations in a sequence.')
parser.add_argument('sequence_file', type=str, help='Path to the input FASTA file')
parser.add_argument('outputdir', type=str, help='Path to the output directory')

# Parse the command-line arguments
args = parser.parse_args()

# Get the sequence file from command-line input
sequence_file = args.sequence_file
outputdir = args.outputdir
# Extract threshold and sdlen from the filename

pattern = r'.*/seq(\d+)_mode\d+_\d+_sim(\d+)_len(\d+)\.fa'

print(sequence_file)
match = re.match(pattern, sequence_file)

if match:
    similarity_threshold = int(match.group(2)) / 100
    sdsize = int(int(match.group(3)))
else:
    raise ValueError("Filename does not match the expected format.")

print(sdsize)
print(similarity_threshold)
# Read the sequence from the file using Biopython
for record in SeqIO.parse(sequence_file, "fasta"):
    sequence = str(record.seq)
    break

# Set parameters
n_mutations = 3  # Number of mutations to introduces

# Get the filename of the input sequence without the extension an without the path

sequence_name = sequence_file.split("/")[1].split(".")[0]



# Apply mutations with checks
mutated_sequence, mutation_log = mutate_sequence(sequence, sdsize, similarity_threshold, n_mutations, sequence_name, outputdir)

# Output the mutated sequence and the log of mutations
#print(mutated_sequence)
for log in mutation_log:
    print(log)