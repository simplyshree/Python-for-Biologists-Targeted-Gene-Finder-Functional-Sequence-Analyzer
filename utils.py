from Bio.Seq import Seq

# GC Content
def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return ((g + c) / len(seq)) * 100

# Transcription (DNA -> RNA)
def transcribe(seq):
    return Seq(seq).transcribe()

# Translation (RNA -> Protein)
def translate(seq):
    seq = seq[:len(seq) - (len(seq) % 3)]  # trim extra bases
    return Seq(seq).translate()

from Bio.Align import PairwiseAligner

def align_sequences(seq1, seq2):
    aligner = PairwiseAligner()
    
    # Set alignment mode
    aligner.mode = 'global'   # like globalxx
    
    alignments = aligner.align(seq1, seq2)
    
    return alignments[0]

# Mutation Detection
def detect_mutations(seq1, seq2):
    mutations = []
    min_len = min(len(seq1), len(seq2))

    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mutations.append((i, seq1[i], seq2[i]))

    if len(seq1) > len(seq2):
        mutations.append(("Insertion", seq1[len(seq2):]))
    elif len(seq2) > len(seq1):
        mutations.append(("Deletion", seq2[len(seq1):]))

    return mutations

# ORF Finder
def find_orfs(seq):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]

            if codon == start_codon:
                for j in range(i, len(seq)-2, 3):
                    stop = seq[j:j+3]
                    if stop in stop_codons:
                        orfs.append((i, j+3))
                        break
            i += 3

    return orfs

print("completed")