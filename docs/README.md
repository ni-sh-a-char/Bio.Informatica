# Bio.Informatica  
**A Common Platform for the study related to Bioinformatics, Computational Biology and Computational Chemistry.**  

---  

## Table of Contents  

| Section | Description |
|---------|-------------|
| **[Overview](#overview)** | What the project is and why it exists |
| **[Installation](#installation)** | How to get the package up and running on your system |
| **[Quick‑Start / Usage](#usage)** | Minimal code snippets to start using the library |
| **[API Documentation](#api-documentation)** | Detailed reference for the public classes, functions and modules |
| **[Examples](#examples)** | Real‑world workflows (sequence analysis, molecular docking, data visualisation, …) |
| **[Contributing](#contributing)** | How to help improve the project |
| **[License & Citation](#license--citation)** | Legal information and how to cite the software |

---  

## Overview  

Bio.Informatica is a **modular, extensible Python library** that brings together common tools needed for modern life‑science research:

* **Bioinformatics** – sequence handling, alignment, annotation, phylogenetics, and NGS data preprocessing.  
* **Computational Biology** – network analysis, systems‑biology modelling, and simulation utilities.  
* **Computational Chemistry** – cheminformatics, molecular descriptor calculation, and simple docking pipelines.  

The library is deliberately **framework‑agnostic**: you can use a single function call for a quick task, or build complex pipelines by chaining the provided components. All heavy‑lifting is delegated to well‑tested third‑party libraries (Biopython, RDKit, NetworkX, scikit‑learn, etc.) while Bio.Informatica supplies the glue code, data‑validation, and a consistent API.

---  

## Installation  

### Prerequisites  

| Tool | Minimum version |
|------|-----------------|
| Python | **3.9** (3.10‑3.12 are fully tested) |
| pip / conda | latest |
| git | optional (for installing from source) |
| C++ compiler | required only for optional RDKit builds (see below) |

> **Tip:** If you plan to use the cheminformatics sub‑module, install **RDKit** first (see the *Optional dependencies* section).

### 1️⃣ Install from PyPI (recommended)

```bash
# Create an isolated environment (highly recommended)
python -m venv bioinf-env
source bioinf-env/bin/activate   # on Windows: .\bioinf-env\Scripts\activate

# Install the core package
pip install bio-informatica
```

### 2️⃣ Install from Conda (if you prefer conda)

```bash
conda create -n bioinf-env python=3.11
conda activate bioinf-env

# The conda‑forge channel hosts the package
conda install -c conda-forge bio-informatica
```

### 3️⃣ Install from source (latest development version)

```bash
# Clone the repository
git clone https://github.com/your-org/Bio.Informatica.git
cd Bio.Informatica

# Install in editable mode with optional extras
pip install -e .[all]   # installs core + all optional dependencies
```

### Optional dependencies  

| Feature | Extra name | Packages installed |
|---------|------------|--------------------|
| **Cheminformatics** | `chem` | `rdkit`, `openbabel` |
| **Machine‑learning pipelines** | `ml` | `scikit-learn`, `xgboost`, `lightgbm` |
| **Visualization** | `viz` | `matplotlib`, `seaborn`, `plotly`, `nglview` |
| **All of the above** | `all` | installs every optional extra |

You can install a specific extra, e.g.:

```bash
pip install bio-informatica[chem]
```

---  

## Usage  

Below is a **minimal “hello‑world”** for each major sub‑module. All examples assume you have activated the environment where Bio.Informatica is installed.

### 1️⃣ Bioinformatics – FASTA parsing & alignment  

```python
from bioinformatica.seq import FastaReader, PairwiseAligner

# Load two sequences from a FASTA file
seqs = FastaReader("data/example.fasta").read()
seq_a, seq_b = seqs["seqA"], seqs["seqB"]

# Perform a global Needleman–Wunsch alignment
aligner = PairwiseAligner(mode="global", match=1, mismatch=-1, gap_open=-2, gap_extend=-0.5)
alignment = aligner.align(seq_a, seq_b)

print(alignment.format())
```

### 2️⃣ Computational Biology – Gene‑regulatory network analysis  

```python
from bioinformatica.network import NetworkBuilder, NetworkAnalyzer

# Build a directed graph from an edge‑list CSV
net = NetworkBuilder.from_edge_table("data/GRN_edges.csv", directed=True)

# Compute centralities and export a summary table
analyzer = NetworkAnalyzer(net)
summary = analyzer.centrality_summary()
summary.to_csv("results/centralities.csv", index=False)
```

### 3️⃣ Computational Chemistry – Molecular descriptor calculation  

```python
from bioinformatica.chem import Molecule, DescriptorCalculator

# Load a SMILES string (or a SDF file)
mol = Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin

# Compute a set of 2D/3D descriptors
calc = DescriptorCalculator()
descriptors = calc.calculate(mol, descriptors=["MolLogP", "TPSA", "NumRotatableBonds"])

print(descriptors)
```

### 4️⃣ End‑to‑end pipeline (FASTA → alignment → phylogeny)  

```python
from bioinformatica.seq import FastaReader, MultipleSeqAligner, PhyloTreeBuilder

# 1. Load a multi‑FASTA file
seqs = FastaReader("data/multiseq.fasta").read_all()

# 2. Perform a progressive multiple sequence alignment
msa = MultipleSeqAligner(method="mafft", options="-auto").align(seqs)

# 3. Build a neighbor‑joining tree
tree = PhyloTreeBuilder(method="nj").build(msa)

# 4. Visualise (requires the `viz` extra)
tree.plot(show=True, output="results/tree.png")
```

---  

## API Documentation  

> The full API reference is also generated automatically by **Sphinx** and hosted at `https://your-org.github.io/Bio.Informatica/`. Below is a concise overview of the most important public classes and functions.

### `bioinformatica.seq`  

| Class / Function | Purpose | Key Parameters |
|------------------|---------|----------------|
| `FastaReader(path: str, encoding: str = "utf-8")` | Parse FASTA files (single or multi‑record). | `path`, `encoding` |
| `SeqRecord` | Light‑weight container for a sequence and its metadata. | `id`, `description`, `seq` |
| `PairwiseAligner(mode: str = "global", **kwargs)` | Wrapper around Biopython’s `PairwiseAligner`. | `mode`, `match`, `mismatch`, `gap_open`, `gap_extend` |
| `MultipleSeqAligner(method: str = "mafft", options: str = "")` | Run external MSA tools (MAFFT, ClustalΩ, MUSCLE). | `method`, `options` |
| `PhyloTreeBuilder(method: str = "nj")` | Build phylogenetic trees from an alignment. | `method`, `bootstrap` |
| `SeqIOHelper` | Convenience functions for format conversion (FASTA ↔ FASTQ ↔ GenBank). | – |

### `bioinformatica.network`  

| Class / Function | Purpose | Key Parameters |
|------------------|---------|----------------|
| `NetworkBuilder` | Construct `networkx.Graph`/`DiGraph` from edge tables, adjacency matrices, or GML files. | `directed`, `weight_col` |
| `NetworkAnalyzer` | Compute topological metrics (degree, betweenness, closeness, eigenvector, clustering). | `graph` |
| `CommunityDetector` | Run community detection algorithms (Louvain, Leiden, Infomap). | `algorithm`, `resolution` |
| `NetworkExporter` | Export graphs to PNG, GraphML, Cytoscape JSON, etc. | `format`, `path` |

### `bioinformatica.chem`  

| Class / Function | Purpose | Key Parameters |
|------------------|---------|----------------|
| `Molecule` | Thin wrapper around RDKit’s `Mol` object with extra I/O helpers. | `from_smiles`, `from_sdf`, `from_mol2` |
| `DescriptorCalculator` | Compute 1‑D, 2‑D, 3‑D descriptors and fingerprints. | `descriptors`, `fingerprints` |
| `DockingRunner` | Simple wrapper for AutoDock Vina (requires Vina binary). | `receptor`, `ligand`, `center`, `size`, `exhaustiveness` |
| `ConformerGenerator` | Generate low‑energy conformers using ETKDG. | `num_confs