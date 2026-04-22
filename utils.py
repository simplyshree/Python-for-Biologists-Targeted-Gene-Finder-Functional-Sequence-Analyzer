from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt


def gc_content(seq):
    return ((seq.count("G") + seq.count("C")) / len(seq)) * 100 if len(seq) > 0 else 0



def transcribe(seq):
    return Seq(seq).transcribe()



def translate(seq):
    seq = seq[:len(seq) - (len(seq) % 3)]
    return Seq(seq).translate()



def align_sequences(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(seq1, seq2)
    return alignments[0]



def detect_mutations(seq1, seq2):
    mutations = []
    min_len = min(len(seq1), len(seq2))

    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mutations.append((i, seq1[i], seq2[i]))

    return mutations



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



def filter_orfs(orfs, min_length=5):
    return [orf for orf in orfs if (orf[1] - orf[0]) >= min_length]



def translate_orfs(seq, orfs):
    proteins = []
    for start, end in orfs:
        coding_seq = seq[start:end]
        coding_seq = coding_seq[:len(coding_seq) - (len(coding_seq) % 3)]
        proteins.append((start, end, Seq(coding_seq).translate()))
    return proteins



def score_orf(seq, start, end):
    length = end - start
    sub_seq = seq[start:end]
    gc = (sub_seq.count("G") + sub_seq.count("C")) / length
    return length * gc


# Mutations inside ORFs
def mutations_in_orfs(mutations, orfs):
    result = []
    for pos, a, b in mutations:
        for start, end in orfs:
            if start <= pos < end:
                result.append((pos, a, b))
    return result



def classify_mutation(seq1, seq2, pos):
    codon_start = (pos // 3) * 3

    codon1 = seq1[codon_start:codon_start+3]
    codon2 = seq2[codon_start:codon_start+3]

    if len(codon1) < 3 or len(codon2) < 3:
        return "Invalid"

    aa1 = Seq(codon1).translate()
    aa2 = Seq(codon2).translate()

    if aa1 == aa2:
        return "Silent"
    elif aa2 == "*":
        return "Nonsense"
    else:
        return "Missense"



def plot_gc_content(seq):
    if len(seq) == 0:
        print("Empty sequence")
        return

    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
    at = 100 - gc

    plt.figure()
    plt.bar(['GC Content', 'AT Content'], [gc, at])
    plt.title("GC vs AT Content")
    plt.ylabel("Percentage")
    plt.show()



def plot_orf_lengths(orfs):
    lengths = [end - start for start, end in orfs]

    if not lengths:
        print("No ORFs to plot")
        return

    plt.figure()
    plt.hist(lengths)
    plt.title("ORF Length Distribution")
    plt.xlabel("Length (bp)")
    plt.ylabel("Frequency")
    plt.show()


print("completed")