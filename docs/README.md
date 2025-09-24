# Bio.Informatica  
**A Common Platform for the study related to Bioinformatics, Computational Biology and Computational Chemistry.**  

---  

## Table of Contents  

| Section | Description |
|---------|-------------|
| **[Installation](#installation)** | How to get Bio.Informatica up and running on your machine (pip, conda, Docker, source). |
| **[Quick‑Start Usage](#quick-start-usage)** | Command‑line interface (CLI) and Python API basics. |
| **[API Documentation](#api-documentation)** | Overview of the public modules, classes, and functions (auto‑generated with Sphinx). |
| **[Examples & Tutorials](#examples--tutorials)** | Ready‑to‑run notebooks and scripts for common bio‑/computational‑chemistry tasks. |
| **[Contributing & Development](#contributing--development)** | Guidelines for extending the platform. |
| **[License & Citation](#license--citation)** | Legal and citation information. |

---  

## Installation  

Bio.Informatica is distributed as a pure‑Python package with optional compiled extensions for heavy‑weight tasks (e.g., molecular dynamics, GPU‑accelerated alignment). Choose the installation method that best fits your workflow.

### 1. Prerequisites  

| Tool | Minimum version | Why it’s needed |
|------|----------------|-----------------|
| Python | **3.9** (≥3.9, ≤3.12) | Core language |
| pip | **23.0** | Package manager |
| conda (optional) | **4.12** | Binary dependencies (e.g., BLAST, OpenMM) |
| Docker (optional) | **20.10** | Containerised reproducibility |
| Git | **2.30** | Source checkout |

> **Tip:** If you plan to use the GPU‑accelerated modules, make sure you have a compatible CUDA toolkit (≥11.8) and the appropriate drivers installed.

### 2. Installing via `pip` (recommended)  

```bash
# Create an isolated environment (highly recommended)
python -m venv bioinf-env
source bioinf-env/bin/activate   # on Windows: .\bioinf-env\Scripts\activate

# Upgrade pip and install the package
pip install --upgrade pip
pip install bioinformatica
```

The default `pip` install pulls **only the pure‑Python core**. To enable optional scientific back‑ends:

```bash
# Install with all optional dependencies (BLAST, OpenMM, RDKit, etc.)
pip install "bioinformatica[all]"
```

### 3. Installing via `conda`  

```bash
conda create -n bioinf-env python=3.11
conda activate bioinf-env

# Bio.Informatica is available on conda-forge
conda install -c conda-forge bioinformatica
```

Conda will automatically resolve binary dependencies (e.g., `openmm`, `rdkit`, `blast`).

### 4. Installing from source (development mode)  

```bash
# Clone the repository
git clone https://github.com/yourorg/Bio.Informatica.git
cd Bio.Informatica

# Install in editable mode with all extras
pip install -e ".[dev,all]"
```

> **Note:** The `-e` flag links the source directory to your environment, so any code changes are reflected instantly.

### 5. Docker image (for reproducible pipelines)  

```bash
docker pull yourorg/bioinformatica:latest
docker run -it --rm \
    -v $(pwd):/workdir \
    -w /workdir \
    yourorg/bioinformatica:latest \
    bash
```

The Docker image contains **Python 3.11**, all optional C‑extensions, and a pre‑installed suite of external tools (BLAST+, GROMACS, OpenMM, RDKit, etc.).

---  

## Quick‑Start Usage  

Bio.Informatica ships both a **command‑line interface (CLI)** and a **Python API**. Below are the most common entry points.

### 1. CLI Overview  

```bash
bioinf --help
```

```
Usage: bioinf [OPTIONS] COMMAND [ARGS]...

  Bio.Informatica – a unified platform for bio‑/computational‑chemistry.

Options:
  -v, --verbose          Increase output verbosity.
  --config FILE          Path to a custom YAML configuration file.
  --version              Show the version and exit.
  --help                 Show this message and exit.

Commands:
  align          Run sequence/structure alignment.
  dock           Perform ligand‑receptor docking (OpenMM/RDKit back‑end).
  mdrun          Run molecular dynamics simulations.
  plot           Quick visualisation of trajectories, spectra, etc.
  workflow       Execute a saved workflow (YAML/JSON).
```

#### Example: Running a BLAST search  

```bash
bioinf align blast \
    --query examples/seqs/fasta1.fasta \
    --db   /data/ncbi/nt \
    --out  results/blast_output.tsv \
    --evalue 1e-5 \
    --max-target-seqs 10
```

#### Example: Docking a small molecule  

```bash
bioinf dock rdkit \
    --receptor examples/pdbs/1abc.pdb \
    --ligand  examples/sdf/ligand.sdf \
    --out     results/docking.sdf \
    --exhaustiveness 12
```

### 2. Python API (core package)  

```python
>>> import bioinformatica as bio
>>> bio.__version__
'2.4.1'
```

#### a. Sequence utilities  

```python
from bio.seq import FastaReader, blast_search

# Load sequences
seqs = FastaReader("examples/seqs/fasta1.fasta").read()

# Run a local BLAST (requires BLAST+ installed)
hits = blast_search(
    query=seqs[0],
    db="/data/ncbi/nt",
    evalue=1e-5,
    max_target_seqs=5,
)
print(hits.head())
```

#### b. Molecular‑dynamics helper  

```python
from bio.md import MDRunner, SystemBuilder

# Build a solvated system from a PDB file
system = SystemBuilder.from_pdb("examples/pdbs/1abc.pdb") \
                      .add_solvent(buffer=10.0) \
                      .add_ions('Na+', 0.15) \
                      .build()

# Run a 10‑ns production simulation (OpenMM back‑end)
runner = MDRunner(system, platform="CUDA")
traj = runner.run(steps=5_000_000, timestep=2.0)   # 10 ns @ 2 fs

traj.save("results/1abc_10ns.dcd")
```

#### c. Chemistry utilities (RDKit wrapper)  

```python
from bio.chem import Molecule, DockingEngine

mol = Molecule.from_sdf("examples/sdf/ligand.sdf")
receptor = Molecule.from_pdb("examples/pdbs/1abc.pdb")

engine = DockingEngine(method="vina")
poses = engine.dock(ligand=mol, receptor=receptor, exhaustiveness=8)

poses[0].to_sdf("results/best_pose.sdf")
```

---  

## API Documentation  

The API is generated automatically with **Sphinx** and hosted on ReadTheDocs (https://bioinformatica.readthedocs.io). Below is a high‑level overview of the most important modules.  

> **Tip:** Use `help(bio.<module>)` or `pydoc` for interactive exploration.

| Module | Purpose | Key Classes / Functions |
|--------|---------|--------------------------|
| `bio.seq` | Sequence I/O, alignment, BLAST, HMMER | `FastaReader`, `FastqReader`, `blast_search()`, `hmmer_search()` |
| `bio.align` | Pairwise & multiple alignment (Clustal‑Omega, MAFFT) | `ClustalAligner`, `MafftAligner`, `AlignmentResult` |
| `bio.struct` | 3‑D structure parsing, manipulation, RMSD calculations | `PDBParser`, `Structure`, `rmsd()` |
| `bio.chem` | Chemistry utilities (RDKit, Open Babel) | `Molecule`, `MolGraph`, `DockingEngine` |
| `bio.md` | Molecular‑dynamics workflow (OpenMM, GROMACS wrappers) | `SystemBuilder`, `MDRunner`, `Trajectory` |
| `bio.plot` | Quick visualisation (matplotlib, plotly) | `plot_rmsd()`, `plot_energy()`, `interactive_viewer()` |
| `bio.workflow` | YAML/JSON based pipeline orchestration | `Workflow`, `Step`, `run_workflow()` |
| `bio.utils` | Helper functions (logging, config handling, file utils) | `