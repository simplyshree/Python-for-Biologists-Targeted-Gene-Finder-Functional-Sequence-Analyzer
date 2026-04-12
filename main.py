from utils import *

def main():
    print("DNA Sequence Analyzer")

    seq1 = input("Enter DNA sequence: ").upper()
    seq2 = input("Enter mutated sequence: ").upper()

    # Basic Analysis
    print("\n--- Basic Analysis ---")
    print("Length:", len(seq1))
    print("GC Content:", gc_content(seq1))

    rna = transcribe(seq1)
    protein = translate(rna)

    print("RNA:", rna)
    print("Protein:", protein)

    # Mutation Detection
    print("\n--- Mutation Detection ---")
    mutations = detect_mutations(seq1, seq2)

    if mutations:
        for m in mutations:
            print("Mutation:", m)
    else:
        print("No mutations found")

    # Alignment
    print("\n--- Alignment ---")
    alignment = align_sequences(seq1, seq2)
    print(str(alignment))

    # ORF Finder
    print("\n--- ORF Finder ---")
    orfs = find_orfs(seq1)

    for i, (start, end) in enumerate(orfs):
        print(f"ORF {i+1}: Start={start}, End={end}, Length={end-start}")

if __name__ == "__main__":
    main()