import tkinter as tk
from tkinter import scrolledtext, messagebox
from utils import *

def run_analysis():
    seq1 = entry_seq1.get().upper()
    seq2 = entry_seq2.get().upper()

    if not seq1 or not seq2:
        messagebox.showerror("Error", "Please enter both sequences")
        return

    output.delete(1.0, tk.END)

    
    output.insert(tk.END, "=== Basic Analysis ===\n")
    output.insert(tk.END, f"Length: {len(seq1)}\n")
    output.insert(tk.END, f"GC Content: {gc_content(seq1):.2f}%\n")

    rna = transcribe(seq1)
    output.insert(tk.END, f"RNA: {rna}\n\n")

    
    output.insert(tk.END, "=== Alignment ===\n")
    alignment = align_sequences(seq1, seq2)
    output.insert(tk.END, str(alignment) + "\n")
    output.insert(tk.END, f"Score: {alignment.score}\n\n")

    
    output.insert(tk.END, "=== ORF Detection ===\n")
    orfs = find_orfs(seq1)
    filtered_orfs = filter_orfs(orfs)

    if not filtered_orfs:
        output.insert(tk.END, "No significant ORFs found\n")
        return

    scored_orfs = []
    for start, end in filtered_orfs:
        score = score_orf(seq1, start, end)
        scored_orfs.append((start, end, score))

    scored_orfs.sort(key=lambda x: x[2], reverse=True)

    for i, (start, end, score) in enumerate(scored_orfs[:3]):
        output.insert(tk.END, f"ORF {i+1}: Start={start}, End={end}, Score={score:.2f}\n")

    
    output.insert(tk.END, "\n=== Proteins ===\n")
    proteins = translate_orfs(seq1, [(s, e) for s, e, _ in scored_orfs[:3]])

    for i, (start, end, protein) in enumerate(proteins):
        output.insert(tk.END, f"\nORF {i+1} Protein:\n{protein}\n")

    
    output.insert(tk.END, "\n=== Mutation Analysis ===\n")
    mutations = detect_mutations(seq1, seq2)
    orf_mutations = mutations_in_orfs(mutations, filtered_orfs)

    if not orf_mutations:
        output.insert(tk.END, "No mutations inside ORFs\n")
    else:
        for pos, a, b in orf_mutations:
            mtype = classify_mutation(seq1, seq2, pos)
            output.insert(tk.END, f"Position {pos}: {a} → {b} ({mtype})\n")



root = tk.Tk()
root.title("DNA Sequence Analyzer")
root.geometry("750x650")


tk.Label(root, text="Original DNA Sequence").pack()
entry_seq1 = tk.Entry(root, width=90)
entry_seq1.pack()

tk.Label(root, text="Mutated DNA Sequence").pack()
entry_seq2 = tk.Entry(root, width=90)
entry_seq2.pack()


tk.Button(root, text="Analyze", command=run_analysis, bg="blue", fg="white").pack(pady=10)

tk.Button(root, text="Show GC Graph",
          command=lambda: plot_gc_content(entry_seq1.get().upper())).pack(pady=5)

tk.Button(root, text="Show ORF Graph",
          command=lambda: plot_orf_lengths(filter_orfs(find_orfs(entry_seq1.get().upper())))).pack(pady=5)


output = scrolledtext.ScrolledText(root, width=90, height=25)
output.pack()

root.mainloop()