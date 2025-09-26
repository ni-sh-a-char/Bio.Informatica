# Bio.Informatica  
**A Common Platform for the study related to Bioinformatics, Computational Biology and Computational Chemistry.**  

---  

## Table of Contents  

| Section | Description |
|---------|-------------|
| **[Overview](#overview)** | What the project is and why it exists |
| **[Installation](#installation)** | How to get the library up and running (pip, conda, source) |
| **[Quick‑Start / Usage](#usage)** | Minimal code to start using the core functionality |
| **[API Documentation](#api-documentation)** | Detailed reference for the public modules, classes and functions |
| **[Examples](#examples)** | Real‑world notebooks & scripts that showcase typical workflows |
| **[Contributing & Development](#contributing--development)** | How to contribute, run tests, and build the docs |
| **[License & Citation](#license--citation)** | Legal information and how to cite the software |

---  

## Overview  

Bio.Informatica is a **modular, extensible Python platform** that brings together tools for:

* **Sequence analysis** – alignment, motif discovery, variant annotation.  
* **Structural bioinformatics** – protein‑ligand docking, molecular dynamics preprocessing, cheminformatics utilities.  
* **Systems biology** – pathway enrichment, network analysis, simulation of kinetic models.  

All components share a **common data model** (based on `pydantic`/`dataclasses`) and a **unified I/O layer** that can read/write FASTA, FASTQ, PDB, SDF, CSV, JSON, and HDF5 files.  

The library is deliberately **lightweight** (core dependencies < 30 MB) and **plug‑in friendly** – you can add custom analysis modules without touching the core code base.

---  

## Installation  

### 1. From PyPI (recommended)

```bash
# Create a clean environment (optional but recommended)
python -m venv bioinf-env
source bioinf-env/bin/activate   # on Windows: .\bioinf-env\Scripts\activate

# Install the latest stable release
pip install bio-informatica
```

> **Tip:** The package name on PyPI is `bio-informatica` (dash, not dot).  

### 2. From Conda (community channel)

```bash
conda create -n bioinf-env python=3.11
conda activate bioinf-env
conda install -c conda-forge bio-informatica
```

### 3. From source (development version)

```bash
# Clone the repository
git clone https://github.com/yourorg/Bio.Informatica.git
cd Bio.Informatica

# Install in editable mode with all optional extras
pip install -e .[all]

# Or install only the core plus selected extras
# pip install -e .[seq,chem]   # seq = sequence tools, chem = cheminformatics
```

#### Optional dependencies (extras)

| Extra | What it adds | Install command |
|-------|--------------|-----------------|
| `seq` | Biopython, pysam, scikit‑bio | `pip install bio-informatica[seq]` |
| `chem`| RDKit, Open Babel, pybel | `pip install bio-informatica[chem]` |
| `ml`  | scikit‑learn, pytorch, tensorflow | `pip install bio-informatica[ml]` |
| `all` | Everything above | `pip install bio-informatica[all]` |

### 4. System requirements  

| Requirement | Minimum version |
|-------------|-----------------|
| Python      | 3.9+ (3.11 recommended) |
| gcc / clang | 7.0+ (for compiled extensions) |
| libstdc++   | 6.0+ |
| CUDA (optional, for ML modules) | 11.2+ |

---  

## Usage  

Below is a **minimal “Hello‑World”** that demonstrates the three pillars of the platform:

```python
# -------------------------------------------------
# 1️⃣  Import the high‑level façade
# -------------------------------------------------
from bioinformatica import Sequence, Molecule, Pathway

# -------------------------------------------------
# 2️⃣  Load a FASTA file and compute a simple statistic
# -------------------------------------------------
seq = Sequence.from_fasta("data/example.fasta")
print(f"Sequence length: {len(seq)}")
print(f"GC content: {seq.gc_content():.2%}")

# -------------------------------------------------
# 3️⃣  Load a small molecule (SDF) and compute a descriptor
# -------------------------------------------------
mol = Molecule.from_sdf("data/ligand.sdf")
print(f"Molecular weight: {mol.molecular_weight():.2f} Da")
print(f"LogP (XLogP3): {mol.logp():.2f}")

# -------------------------------------------------
# 4️⃣  Perform a quick pathway enrichment
# -------------------------------------------------
genes = ["TP53", "BRCA1", "EGFR", "MYC"]
enrichment = Pathway.enrich(genes, db="KEGG")
print(enrichment.head())
```

### Command‑line interface (CLI)

Bio.Informatica ships with a small but handy CLI (`bioinf`).  

```bash
# Show help
bioinf --help

# Convert a FASTQ to FASTA
bioinf seq convert --in reads.fastq --out reads.fasta

# Compute basic molecular descriptors for a whole SDF library
bioinf chem descriptors --in library.sdf --out descriptors.csv

# Run a quick KEGG enrichment on a gene list
bioinf pathway enrich --genes genes.txt --db kegg --out enrichment.tsv
```

All CLI commands accept a `--verbose` flag for detailed logging and a `--threads N` flag to parallelise where possible.

---  

## API Documentation  

> **Note:** The documentation below reflects the **v2.3.0** release. Use `help(bioinformatica)` or `pydoc` for the most up‑to‑date signatures.

### Core Packages  

| Package | Description | Most important objects |
|---------|-------------|------------------------|
| `bioinformatica.sequence` | Utilities for nucleic‑acid and protein sequences. | `Sequence`, `Alignment`, `MotifFinder` |
| `bioinformatica.chemistry` | Cheminformatics, molecular descriptors, file conversion. | `Molecule`, `DockingEngine`, `Fingerprint` |
| `bioinformatica.pathway` | Enrichment, network analysis, simulation. | `Pathway`, `Network`, `FluxSimulator` |
| `bioinformatica.io` | Unified I/O layer (auto‑detects format). | `read`, `write`, `FileHandler` |
| `bioinformatica.ml` | Machine‑learning wrappers (classification, clustering). | `Classifier`, `Regressor`, `Embedding` |
| `bioinformatica.utils` | Helper functions (logging, tqdm wrappers, caching). | `logger`, `cached_property` |

---

### 1. `bioinformatica.sequence`

```python
class Sequence:
    """Container for a single biological sequence (DNA, RNA or protein)."""

    @classmethod
    def from_fasta(cls, path: str, *, id: str | None = None) -> "Sequence": ...
    @classmethod
    def from_fastq(cls, path: str) -> "Sequence": ...

    def __len__(self) -> int: ...
    def gc_content(self) -> float: ...
    def translate(self, *, to_stop: bool = True) -> "Sequence": ...

    def reverse_complement(self) -> "Sequence": ...
    def subseq(self, start: int, end: int) -> "Sequence": ...

    # Alignment helpers
    def align(self, other: "Sequence", *, method: str = "global") -> "Alignment": ...
```

```python
class Alignment:
    """Result of a pairwise or multiple alignment."""

    score: float
    aligned_seq1: str
    aligned_seq2: str
    start1: int
    start2: int

    def identity(self) -> float: ...
    def to_fasta(self, path: str) -> None: ...
```

---

### 2. `bioinformatica.chemistry`

```python
class Molecule:
    """Thin wrapper around RDKit's Mol object with extra utilities."""

    @classmethod
    def from_sdf(cls, path: str, *, index: int = 0) -> "Molecule": ...
    @classmethod
    def from_smiles(cls, smiles: str) -> "Molecule": ...

    def molecular_weight(self) -> float: ...
    def logp(self) -> float: ...
    def tpsa(self) -> float: ...
    def fingerprint(self, *, radius: int = 2, nbits: int = 2048) -> "Fingerprint": ...

    # Simple docking interface (requires AutoDock Vina or OpenBabel)
    def dock(self, receptor: "Molecule", *, exhaustiveness: int = 8) -> "DockingResult": ...
```

```python
class Fingerprint:
    bits: bytes
    nbits: int

    def tanimoto(self, other: "Fingerprint") -> float: ...
    def to_hex(self) -> str: ...
```

---

### 3. `bioinformatica.pathway`

```python