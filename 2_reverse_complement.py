from Bio.Seq import Seq
from Bio import SeqIO

# Function to calculate reverse complement
def reverse_complement(sequence):
    reverse_comp = sequence.reverse_complement()
    return reverse_comp

# Main function to handle user's choice of input
def main():
    print("Choose input method:")
    print("1. Manual sequence input")
    print("2. FASTA file input")
    choice = input("Enter your choice (1 or 2): ")

    # Manual sequence input
    if choice == '1':
        dna_sequence = Seq(input("Enter the DNA sequence: "))
        reverse_comp_sequence = reverse_complement(dna_sequence)
        print(f"Reverse Complement: {reverse_comp_sequence}")

    # FASTA file input
    elif choice == '2':
        fasta_file = input("Enter the path to your FASTA file: ")
        try:
            # Read the FASTA file using SeqIO
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            for i, record in enumerate(sequences):
                # Calculate the reverse complement for each sequence in the FASTA file
                reverse_comp_sequence = reverse_complement(record.seq)
                print(f"Reverse Complement for Sequence {i+1} ({record.id}): {reverse_comp_sequence}")
        except FileNotFoundError:
            print("Error: The specified FASTA file was not found.")
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        print("Invalid choice. Please select 1 or 2.")

# Call the main function
if __name__ == "__main__":
    main()
