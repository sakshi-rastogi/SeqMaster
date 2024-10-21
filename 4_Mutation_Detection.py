from Bio import SeqIO

def detect_mutations(seq1, seq2):
    mutations = []
    length = max(len(seq1), len(seq2))

    for i in range(length):
        # Check for substitutions
        if i < len(seq1) and i < len(seq2):
            if seq1[i] != seq2[i]:
                mutations.append((i + 1, 'substitution', seq1[i], seq2[i]))
        
        # Check for insertions or deletions
        elif i < len(seq1):  # Deletion in seq2
            mutations.append((i + 1, 'deletion', seq1[i], '-'))
        elif i < len(seq2):  # Insertion in seq1
            mutations.append((i + 1, 'insertion', '-', seq2[i]))

    return mutations

def main():
    # User choice for input
    choice = input("Do you want to enter sequences manually or use FASTA files? (manual/fasta): ").strip().lower()

    if choice == 'manual':
        seq1 = input("Enter the first DNA sequence: ")
        seq2 = input("Enter the second DNA sequence: ")
    elif choice == 'fasta':
        fasta_file = input("Enter the path to your FASTA file (two sequences): ")
        sequences = list(SeqIO.parse(fasta_file, "fasta"))

        # Assuming the first two sequences in the FASTA file
        seq1 = str(sequences[0].seq)
        seq2 = str(sequences[1].seq)
    else:
        print("Invalid choice. Exiting.")
        return

    # Detect mutations
    mutations = detect_mutations(seq1, seq2)

    # Print results
    if mutations:
        print("Mutations detected:")
        for mutation in mutations:
            position, mutation_type, original, new = mutation
            print(f"Position: {position}, Type: {mutation_type}, Original: {original}, New: {new}")
    else:
        print("No mutations detected.")

if __name__ == "__main__":
    main()
