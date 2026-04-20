# Targeted Gene Finder & Functional Sequence Analyzer
# 🧬 DNA Sequence Analyzer

### Targeted Gene Prediction and Mutation Impact Analysis

## Overview

This project is a bioinformatics tool developed using Python that performs DNA sequence analysis, including ORF detection, mutation analysis, transcription, translation, and graphical visualization.

The system is built using **Biopython** and includes a user-friendly GUI for easy interaction.

---

## Features

* DNA → RNA transcription
* RNA → Protein translation
* Open Reading Frame (ORF) detection
* Mutation detection and classification
* Sequence alignment
* GC content analysis
* Graph visualization (GC content & ORF distribution)
* GUI-based interface

---

## Technologies Used

* Python
* Biopython
* Tkinter
* Matplotlib

---

## 📂 Project Structure

```
DNA_Analyzer/
│── utils.py     # Core logic
│── main.py      # CLI version
│── gui.py       # GUI version
```

---

## ▶How to Run

### 1. Install Dependencies

```
pip install biopython matplotlib
```

### 2. Run GUI Version (Recommended)

```
python gui.py
```

### 3. Run CLI Version

```
python main.py
```

---

## Sample Input

### Original DNA:

```
ATGAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCTAA
```

### Mutated DNA:

```
ATGAAAGGGTTTCCCAAAGGGTCTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCAAAGGGTTTCCCTAA
```

---

## Output

* ORF positions and scores
* Protein sequences
* Mutation classification
* Graphical plots

---


---

##  References

* Biopython Documentation
* Python Documentation
* NCBI Resources

---

## 🏁 Conclusion

This project demonstrates how computational tools can be used to analyze biological sequences efficiently and visualize important genetic information.

---

