from Bio import SeqIO

def is_palindrome(seq):
    """Check if a given sequence is a palindrome."""
    return seq == seq[::-1]

def find_palindromes(sequence, min_length=4):
    """Find palindromic sequences in the given DNA sequence."""
    palindromes = []
    seq_length = len(sequence)

    # Check for palindromes of length min_length to seq_length
    for length in range(min_length, seq_length + 1):
        for start in range(seq_length - length + 1):
            subseq = sequence[start:start + length]
            if is_palindrome(subseq):
                palindromes.append(subseq)

    return palindromes

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

    # Find palindromes
    palindromes = find_palindromes(dna_sequence)

    # Print the results
    if palindromes:
        print("Palindromic sequences found:")
        for palindrome in set(palindromes):  # Use set to avoid duplicates
            print(palindrome)
    else:
        print("No palindromic sequences found.")

if __name__ == "__main__":
    main()
