# 🧬 PCR Primer Design & Analysis Pipeline

A small, reproducible bioinformatics workflow that goes from a **GENCODE transcript FASTA file** to a **validated PCR primer pair**: parsing sequences with Biopython, designing primers with NCBI Primer-BLAST, and checking primer quality with the Sequence Manipulation Suite (SMS).

> 🎥 A full screen-recorded walkthrough of the project is available in the demo video: `Primer_Analysis.mp4` *(https://lnkd.in/p/eVXT7tUT)*.

---

## 📌 Table of Contents

- [Overview](#-overview)
- [Workflow](#-workflow)
- [Tools & Data](#-tools--data)
- [Step-by-Step](#-step-by-step)
- [Results](#-results)
- [How to Reproduce](#-how-to-reproduce)
- [Limitations & Future Work](#-limitations--future-work)
- [Author](#-author)

---

## 🔎 Overview

Designing good PCR primers is a core step in molecular biology experiments (gene expression, cloning, validation). A good primer pair should have:

- A suitable **length** (about 18–25 nt)
- **GC content** around 40–60%
- A **melting temperature (Tm)** close to the reaction conditions, with both primers having similar Tm
- Low **self-complementarity** (to avoid primer-dimers and hairpins)
- **Specificity** to the target transcript only

This project demonstrates an end-to-end process to extract a target transcript, design primers for it, and evaluate their quality computationally before ordering them in the lab.

---

## 🔁 Workflow

```
GENCODE FASTA (protein-coding transcripts, v50)
          │
          ▼
Biopython (SeqIO) ── parse & select a target transcript
          │
          ▼
NCBI Primer-BLAST ── design & check specificity of primer pairs
          │
          ▼
SMS: PCR Primer Stats ── Tm, GC%, MW, self-annealing, hairpin checks
          │
          ▼
Final primer pair selection
```

---

## 🛠 Tools & Data

| Category | Details |
|---|---|
| Language | Python 3 (virtual environment `bio_env`) |
| Library | [Biopython](https://biopython.org/) (`Bio.SeqIO`) |
| Environment | Ubuntu (WSL) terminal |
| Input data | `gencode.v50.pc_transcripts.fa` — GENCODE v50 protein-coding transcript sequences |
| Primer design | [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) (Primer3 + BLAST) |
| Primer QC | [Sequence Manipulation Suite – PCR Primer Stats](https://www.bioinformatics.org/sms2/pcr_primer_stats.html) |

---

## 🧭 Step-by-Step

### 1. Parse the GENCODE FASTA file

```python
from Bio import SeqIO

filename = 'gencode.v50.pc_transcripts.fa'
seq = SeqIO.parse(filename, 'fasta')

# Load all records into a list
L = []
for record in seq:
    L.append(record)

# Inspect a transcript
print(L[3].name)
print(L[3].seq)
```

Each GENCODE header packs several fields separated by `|` (transcript ID, gene ID, gene name, transcript length, and UTR/CDS coordinates). For the selected transcript:

| Field | Value |
|---|---|
| Transcript ID | `ENST00000114008.1` |
| Gene ID | `ENSG00000187634.15` |
| Gene name | `SAMD11` (transcript `SAMD11-219`) |
| Transcript length | 2523 nt |
| 5′ UTR | 1–537 |
| CDS | 538–1857 |
| 3′ UTR | 1858–2523 |

### 2. Design primers with NCBI Primer-BLAST

The transcript sequence was pasted as the **PCR template**, with these main parameters:

| Parameter | Value |
|---|---|
| Forward primer search range | 1 – 400 |
| Reverse primer search range | 2000 – 2406 |
| Number of primers to return | 10 |
| Primer Tm (min / opt / max) | 57.0 / 60.0 / 63.0 °C |
| Max Tm difference | 3 °C |

Primer-BLAST returned 10 candidate primer pairs with their position, length, Tm, GC%, self-complementarity and product length.

### 3. Evaluate primer quality with PCR Primer Stats

The forward and reverse primers of the chosen pair were submitted to SMS **PCR Primer Stats** to get detailed properties (base counts, molecular weight, Tm estimates) and a pass/warning report for several PCR suitability tests.

---

## 📊 Results

### Candidate pairs from Primer-BLAST (first five)

| Pair | Forward (5′→3′) | Reverse (5′→3′) | Fwd Tm | Rev Tm | Product (bp) |
|---|---|---|---|---|---|
| 1 | `GGATCGATAGCAGCTCCATGT` | `TGGCTGCATAGACGACTCAC` | 59.73 | 59.83 | 1803 |
| 2 | `TCCGGGATCGATAGCAGCTC` | `CTGGCTGCATAGACGACTCA` | 61.24 | 59.54 | 1808 |
| 3 | `GATAGCAGCTCCATGTCTCCG` | `ACTGGCTGCATAGACGACTC` | 60.34 | 59.54 | 1800 |
| 4 | `CGGGATCGATAGCAGCTCCA` | `CCACTGGCTGCATAGACGAC` | 61.53 | 60.81 | 1809 |
| 5 | `TCGATAGCAGCTCCATGTCTC` | `GTCCACTGGCTGCATAGACG` | 59.39 | 60.81 | 1806 |

*All five pairs fall in the same region, so product sizes are very similar (about 1.8 kb).*

### Detailed QC for Pair 1 (SMS PCR Primer Stats)

| Property | Forward primer | Reverse primer |
|---|---|---|
| Length | 21 nt | 20 nt |
| GC content | 52.38 % | 55.00 % |
| Base counts | G=6, A=5, T=5, C=5 | G=5, A=5, T=4, C=6 |
| Molecular weight | 6446.25 Da | 6102.02 Da |
| Tm (nearest neighbor) | 64.93 °C | 65.09 °C |
| Single base runs / Dinucleotide runs | Pass | Pass |
| Length / %GC / GC clamp | Pass | Pass |
| Hairpin formation | Pass | Pass |
| Tm check (nearest neighbor) | ⚠️ Warning (> 58 °C) | ⚠️ Warning (> 58 °C) |
| Self-annealing | ⚠️ Warning (more than 3 self-annealing bases in a row) | Pass |

### Interpretation

- Both primers have a good length and GC content, no problematic base runs, and no hairpins.
- The **Tm warning** appears because SMS uses a different Tm calculation (nearest-neighbor with its own salt and primer-concentration defaults: 50 mM K⁺/Na⁺, 1.5 mM Mg²⁺, 200 nM primer) than Primer-BLAST, so the values differ between tools. Always match the Tm method to the polymerase/buffer you will use.
- The **forward primer** shows a short self-annealing stretch; this may promote primer-dimer formation. Other pairs from the Primer-BLAST list can be tested as alternatives.

---

## ▶️ How to Reproduce

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>

# 2. Create and activate a virtual environment
python3 -m venv bio_env
source bio_env/bin/activate

# 3. Install dependencies
pip install biopython

# 4. Download the GENCODE protein-coding transcripts FASTA (v50)
#    from https://www.gencodegenes.org/ and place it in the project folder
#    (file name: gencode.v50.pc_transcripts.fa)

# 5. Run the parsing step
python Sequences.py
```

Then continue manually with Primer-BLAST and PCR Primer Stats as described in [Step-by-Step](#-step-by-step).

---

## 🚧 Limitations & Future Work

- Primer-BLAST and SMS steps are done through web tools; they could be automated (for example, with NCBI's web services or the `primer3-py` package).
- Specificity was assessed through Primer-BLAST only; experimental validation (PCR and gel electrophoresis) is still needed.
- Planned improvements:
  - [ ] Automate primer design for any transcript ID or gene name
  - [ ] Compute Tm, GC% and self-complementarity locally in Python
  - [ ] Add a ranking score to pick the best primer pair automatically
  - [ ] Export results to CSV and add plots

---

## 👤 Author

**<Hozaifa Alzohare>**
Biotechnology | Bioinformatics | Ai in Healthcare

---

## 📚 References

- GENCODE — https://www.gencodegenes.org/
- NCBI Primer-BLAST — https://www.ncbi.nlm.nih.gov/tools/primer-blast/
- Sequence Manipulation Suite — https://www.bioinformatics.org/sms2/
- Biopython — https://biopython.org/
