from Bio.Seq import Seq
from Bio import SeqIO

# Function to calculate GC content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100
    return gc_content

# Ask the user how they want to input the sequences
choice = input("Do you want to input sequences manually or from a FASTA file? (manual/fasta): ").lower()

if choice == 'manual':
    # Single or multiple sequence choice
    manual_choice = input("Do you want to analyze a single sequence or multiple sequences? (single/multiple): ").lower()

    if manual_choice == 'single':
        # Input a single sequence
        dna_sequence = Seq(input("Enter DNA sequence: "))
        gc_content = calculate_gc_content(dna_sequence)
        print(f"GC Content: {gc_content:.2f}%")

    elif manual_choice == 'multiple':
        # Ask how many sequences the user wants to analyze
        num_sequences = int(input("How many sequences do you want to analyze? "))
        sequences = []
        for i in range(1, num_sequences + 1):
            seq_input = Seq(input(f"Enter your Sequence {i}: "))
            sequences.append(seq_input)
        
        # Calculate and print GC content for each sequence
        for i, seq in enumerate(sequences):
            gc_content = calculate_gc_content(seq)
            print(f"GC Content for Sequence {i+1}: {gc_content:.2f}%")

elif choice == 'fasta':
    # Ask for the path to the FASTA file
    fasta_file = input("Enter the path to your FASTA file: ")

    # Read the FASTA file and calculate GC content for each sequence
    try:
        # Use Biopython's SeqIO to parse the FASTA file
        sequences = list(SeqIO.parse(fasta_file, "fasta"))

        # Loop through each sequence record in the file
        for i, record in enumerate(sequences):
            # record.seq contains the sequence
            gc_content = calculate_gc_content(record.seq)
            print(f"GC Content for Sequence {i+1} ({record.id}): {gc_content:.2f}%")
    
    except FileNotFoundError:
        print("Error: The specified FASTA file was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

else:
    print("Invalid input choice.")
