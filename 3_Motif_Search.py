from Bio import SeqIO

def search_motif(sequence, motif):
    """Search for the motif in the given sequence."""
    positions = []
    start = 0
    while True:
        # Find the motif in the sequence
        start = sequence.find(motif, start)
        if start == -1:  # No more occurrences found
            break
        positions.append(start)
        start += 1  # Move one position forward to find the next occurrence
    return positions

def main():
    choice = input("Do you want to enter a sequence manually or use a FASTA file? (manual/fasta): ").lower()
    
    if choice == 'manual':
        dna_sequence = input("Enter your DNA sequence: ").upper()
    elif choice == 'fasta':
        fasta_file = input("Enter the path to your FASTA file: ")
        # Read the FASTA file and get the first sequence
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if sequences:
            dna_sequence = str(sequences[0].seq).upper()  # Convert Seq object to string
        else:
            print("No sequences found in the FASTA file.")
            return
    else:
        print("Invalid choice.")
        return

    motif = input("Enter the motif to search for: ").upper()
    
    # Search for the motif
    positions = search_motif(dna_sequence, motif)
    
    if positions:
        print(f"Motif '{motif}' found at positions: {positions}")
    else:
        print(f"Motif '{motif}' not found in the sequence.")

if __name__ == "__main__":
    main()
