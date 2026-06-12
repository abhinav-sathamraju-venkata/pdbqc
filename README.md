# pdbqc

A Python tool for structure quality reporting and cleaning of PDB/mmCIF files. Detects and reports common structural pathologies — alternate locations, insertion codes, microheterogeneity, missing residues, and selenomethionine conversions — in a clean, colored terminal report.

## Install

```bash
pip install pdbqc
```

## Usage

### Python
```python
from pdbqc import clean_structure, format_report

# Basic usage
qc, structure = clean_structure("my_structure.cif")
format_report(qc)

# With CIF file for coverage and missing residue data
qc, structure = clean_structure("my_structure.pdb", cif_path="my_structure.cif")
format_report(qc)

# Keep ligands in cleaned structure
qc, structure = clean_structure("my_structure.pdb", keep_ligands=True)
format_report(qc)
```

### CLI
```bash
pdbqc my_structure.pdb
pdbqc my_structure.pdb --cif my_structure.cif
pdbqc my_structure.pdb --keep-ligands
```

## Output
╔══════════════════════════════════════╗
║ Structure Quality Report             ║
╚══════════════════════════════════════╝
Atoms: 2671
Altlocs detected: none
Insertion codes: none
Microheterogeneity: none
Chain  Missing Residues  Coverage
──────────────────────────────────
A      15                90.5%
B      13                91.8%
C      0                 100.0%
D      0                 100.0%
MSE residues converted: none

## What it checks

- **Alternate locations** — counts atoms with altlocs before filtering; keeps highest occupancy
- **Insertion codes** — detects and reports non-standard residue numbering
- **Microheterogeneity** — flags residues with conflicting identities across altlocs
- **Missing residues** — compares observed CA atoms to SEQRES (requires CIF)
- **Coverage** — per-chain percentage of observed vs expected residues
- **MSE conversion** — selenomethionine residues renamed to MET for standardized downstream use

## Links

- [PyPI](https://pypi.org/project/pdbqc/)
- [GitHub](https://github.com/abhinav-sathamraju-venkata/pdbqc)