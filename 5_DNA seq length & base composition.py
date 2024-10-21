from Bio import SeqIO

def calculate_length_and_composition(sequence):
    """Calculate the length and base composition of the DNA sequence."""
    length = len(sequence)
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')

    # Calculate base composition percentages
    a_percent = (a_count / length) * 100 if length > 0 else 0
    t_percent = (t_count / length) * 100 if length > 0 else 0
    c_percent = (c_count / length) * 100 if length > 0 else 0
    g_percent = (g_count / length) * 100 if length > 0 else 0

    return length, a_percent, t_percent, c_percent, g_percent

def main():
    choice = input("Do you want to enter a sequence manually or use a FASTA file? (manual/fasta): ").lower()
    
    if choice == 'manual':
        dna_sequence = input("Enter your DNA sequence: ").upper()
    elif choice == 'fasta':
        fasta_file = input("Enter the path to your FASTA file: ")
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if sequences:
            dna_sequence = str(sequences[0].seq).upper()  # Convert Seq object to string
        else:
            print("No sequences found in the FASTA file.")
            return
    else:
        print("Invalid choice.")
        return

    # Calculate length and base composition
    length, a_percent, t_percent, c_percent, g_percent = calculate_length_and_composition(dna_sequence)

    # Print the results
    print(f"Length of the DNA sequence: {length} bases")
    print(f"Base Composition:")
    print(f"A: {a_percent:.2f}%")
    print(f"T: {t_percent:.2f}%")
    print(f"C: {c_percent:.2f}%")
    print(f"G: {g_percent:.2f}%")

if __name__ == "__main__":
    main()
