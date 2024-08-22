import argparse
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

def generate_sequences(n_sequences, output_folder, segment_lengths, segment_similarities, seed=None):
    if seed is not None:
        random.seed(seed)
    
    seq_counter = 0
    os.makedirs(output_folder, exist_ok=True)

    for seg_length in segment_lengths:
        for similarity in segment_similarities:
            total_seq_length = seg_length * 10
            
            for seq_num in range(n_sequences):
                seq_counter += 1
                
                 # Step 1: Choose four positions between 20% and 80% of total length, multiples of seg_length
                segment_positions = sorted(random.sample(range(int(0.0 * total_seq_length), int(1 * total_seq_length), seg_length), 4))
                pos_a, pos_b, pos_c, pos_d = segment_positions

                # Step 2: Choose mode and orientations
                segment_order_mode = random.choice([1, 2])
                orientation_flags = [random.choice([0, 1]) for _ in range(2)]

                # Step 3: Generate 10 random sequences of seg_length
                random_sequences = [''.join(random.choices('ACGT', k=seg_length)) for _ in range(10)]

                # Step 4: Replace sequences as per the mode
                if segment_order_mode == 1:
                    random_sequences[pos_c // seg_length] = mutate_sequence(random_sequences[pos_a // seg_length], similarity, orientation_flags[0])
                    random_sequences[pos_d // seg_length] = mutate_sequence(random_sequences[pos_b // seg_length], similarity, orientation_flags[1])
                else:  # segment_order_mode == 2
                    random_sequences[pos_d // seg_length] = mutate_sequence(random_sequences[pos_a // seg_length], similarity, orientation_flags[0])
                    random_sequences[pos_c // seg_length] = mutate_sequence(random_sequences[pos_b // seg_length], similarity, orientation_flags[1])

                starter_sequence = ''.join(random.choices('ACGT', k=seg_length*5))
                end_sequence = ''.join(random.choices('ACGT', k=seg_length*5))

                # Step 5: Merge sequences into one segment
                combined_sequence = ''.join([''.join([starter_sequence, ''.join(random_sequences)]), end_sequence])
                seq_record = SeqRecord(Seq(combined_sequence), id=f"seq{seq_counter}_mode{segment_order_mode}_orient{str(orientation_flags[0])+str(orientation_flags[1])}_sim{int(similarity*100)}_len{seg_length}", description="")

                # Step 6: Save to file
                output_file = os.path.join(output_folder, f"seq{seq_counter}_mode{segment_order_mode}_{str(orientation_flags[0])+str(orientation_flags[1])}_sim{int(similarity*100)}_len{seg_length}.fa")
                SeqIO.write(seq_record, output_file, "fasta")


def mutate_sequence(sequence, similarity, reverse_complement_flag):
    sequence = Seq(sequence)
    if reverse_complement_flag:
        sequence = sequence.reverse_complement()

    num_snps = int(len(sequence) * (1 - similarity))
    sequence = list(str(sequence))

    for _ in range(num_snps):
        snp_pos = random.randint(0, len(sequence) - 1)
        original_base = sequence[snp_pos]
        sequence[snp_pos] = random.choice([base for base in 'ACGT' if base != original_base])

    return ''.join(sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate sequences with specific lengths and similarities.")
    parser.add_argument("n_sequences", type=int, help="Number of sequences per length-similarity combination.")
    parser.add_argument("output_folder", type=str, help="Output folder to save the generated sequences.")
    parser.add_argument("--segment_lengths", type=int, nargs="+", default=[100, 500, 1000, 10000], help="List of segment lengths.")
    parser.add_argument("--segment_similarities", type=float, nargs="+", default=[0.9, 0.95, 0.99], help="List of segment similarities as percentages (e.g., 0.9, 0.95, 0.99).")
    parser.add_argument("--seed", type=int, default=None, help="Seed for random number generator to ensure reproducibility.")

    args = parser.parse_args()

    generate_sequences(args.n_sequences, args.output_folder, args.segment_lengths, args.segment_similarities, args.seed)

   