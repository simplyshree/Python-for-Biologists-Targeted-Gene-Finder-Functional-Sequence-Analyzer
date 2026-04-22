from utils import *

def main():
    print("=== DNA Sequence Analyzer (CLI Version) ===")

    seq1 = input("Enter original DNA sequence: ").upper()
    seq2 = input("Enter mutated DNA sequence: ").upper()

    if not seq1 or not seq2:
        print("Error: Please enter valid sequences")
        return

    
    print("\n--- Basic Analysis ---")
    print("Length:", len(seq1))
    print("GC Content: {:.2f}%".format(gc_content(seq1)))

    rna = transcribe(seq1)
    print("RNA:", rna)

    
    print("\n--- Alignment ---")
    alignment = align_sequences(seq1, seq2)
    print(alignment)
    print("Score:", alignment.score)

    
    print("\n--- ORF Detection ---")
    orfs = find_orfs(seq1)
    filtered_orfs = filter_orfs(orfs)

    if not filtered_orfs:
        print("No significant ORFs found.")
        return

    
    scored_orfs = []
    for start, end in filtered_orfs:
        score = score_orf(seq1, start, end)
        scored_orfs.append((start, end, score))

    scored_orfs.sort(key=lambda x: x[2], reverse=True)

    for i, (start, end, score) in enumerate(scored_orfs[:3]):
        print(f"ORF {i+1}: Start={start}, End={end}, Length={end-start}, Score={score:.2f}")

    
    print("\n--- Protein Sequences ---")
    proteins = translate_orfs(seq1, [(s, e) for s, e, _ in scored_orfs[:3]])

    for i, (start, end, protein) in enumerate(proteins):
        print(f"\nORF {i+1} Protein:")
        print(protein)

    
    print("\n--- Mutation Analysis ---")
    mutations = detect_mutations(seq1, seq2)
    orf_mutations = mutations_in_orfs(mutations, filtered_orfs)

    if not orf_mutations:
        print("No mutations inside ORFs.")
    else:
        for pos, a, b in orf_mutations:
            mtype = classify_mutation(seq1, seq2, pos)
            print(f"Position {pos}: {a} → {b} ({mtype})")


if __name__ == "__main__":
    main()