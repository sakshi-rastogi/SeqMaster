from Bio.Seq import Seq
from Bio.SeqUtils import seq3

# Genetic code dictionary
genetic_code = {
    'AUG': 'M',  # Start codon
    'UUU': 'F', 'UUC': 'F',  # Phenylalanine
    'UUA': 'L', 'UUG': 'L',  # Leucine
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',  # Serine
    'UAU': 'Y', 'UAC': 'Y',  # Tyrosine
    'UGU': 'C', 'UGC': 'C',  # Cysteine
    'UGG': 'W',  # Tryptophan
    'CUU': 'L', 'CUC': 'L', 'CCA': 'P', 'CCU': 'P', 'CCG': 'P',  # Proline
    'AAC': 'N', 'AAU': 'N',  # Asparagine
    'CCA': 'Q', 'CAG': 'Q',  # Glutamine
    'GAA': 'E', 'GAG': 'E',  # Glutamic acid
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',  # Isoleucine
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',  # Valine
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',  # Stop codons
}

def translate_dna_to_protein(dna_sequence):
    # Transcribe DNA to mRNA
    mRNA_sequence = dna_sequence.replace('T', 'U')

    # Translate mRNA to protein
    protein_sequence = ""
    for i in range(0, len(mRNA_sequence), 3):  # Process each codon (3 nucleotides)
        codon = mRNA_sequence[i:i+3]  # Get the codon
        if len(codon) == 3:  # Ensure the codon is complete
            amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codon
            if amino_acid == 'Stop':
                break  # Stop translation at stop codon
            protein_sequence += amino_acid
    
    return protein_sequence

def main():
    # User choice for input
    choice = input("Do you want to enter a DNA sequence manually or use a FASTA file? (manual/fasta): ").strip().lower()

    if choice == 'manual':
        dna_sequence = input("Enter the DNA sequence: ").strip().upper()  # Convert to uppercase for consistency
    elif choice == 'fasta':
        fasta_file = input("Enter the path to your FASTA file: ")
        from Bio import SeqIO
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        dna_sequence = str(sequences[0].seq)  # Assuming the first sequence in the FASTA file
    else:
        print("Invalid choice. Exiting.")
        return

    # Translate DNA to protein
    protein_sequence = translate_dna_to_protein(dna_sequence)

    # Print the protein sequence
    print(f"Protein sequence: {protein_sequence}")

if __name__ == "__main__":
    main()
