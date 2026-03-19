# Interaction Map Explorer (IME)

<p align="center">
  <img src="images/logo_IME.png" alt="IME Logo" />
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.19100311"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19100311.svg" alt="DOI"></a>
</p>

A Python tool that automates the identification and characterization of **Critical Points (CPs)** at the interface between molecular fragments, enabling the study of **intermolecular interactions** in molecule–molecule and protein–ligand systems. It supports Gaussian (`.fchk`) and ORCA (`.wfn`/`.wfx`) wavefunction formats.

Designed for computational chemists who need to systematically map and quantify how fragments interact — whether through hydrogen bonds, van der Waals contacts, or any other type of interaction revealed by topological analysis (QTAIM).

---

## Key Features

- **Automatic CP Detection**: Finds Bond CPs (3,−1), Ring CPs (3,+1) and Cage CPs (3,+3) at the interface between user-defined fragments.
- **Two Execution Modes**:
  - **Interactive Mode**: 3D graphical selector (powered by [FragmentFinder](https://github.com/1JELC1/FragmentFinder)) to manually assign atoms to fragments for each file.
  - **Batch Mode**: Provide one reference ligand `.xyz` file and the tool automatically identifies it in all complexes using graph isomorphism, enabling fully automated batch processing.
- **VMD-Ready Output**: Generates `.tcl` VMD scripts to directly visualize the interface CPs and bond paths on the molecular structure.
- **CP Property Extraction**: Extracts 25 topological properties at each interface CP (density, Laplacian, ELF, LOL, RDG, ESP, etc.) and exports all values to a per-molecule `.csv`.
- **Skip-If-Exists Logic**: Detects already-processed files and skips redundant calculations, allowing safe resumption after interruptions.
- **Flexible Input**: Supports `.fchk`, `.wfn`, and `.wfx` wavefunction files.

---

## Workflow & Methodology

The script follows a structured pipeline:

1. **Directory Configuration**: User provides the wavefunction input folder and output folder paths.
2. **Mode Selection**: User chooses Interactive or Batch mode (ligand-environment).
3. **PDB Generation**: Generates `.pdb` files for the molecular structure, CPs, and bond paths for each wavefunction file.
4. **Fragment Definition**: Fragments can be defined in two ways:
   - *Interactive Mode*: the user manually selects and assigns atoms to different fragments through the 3D graphical interface.
   - *Batch Mode*: the user provides a reference `.xyz` file that defines Fragment 1; all remaining atoms in each complex are automatically assigned to Fragment 2 (environment). The tool then searches for interactions between the ligand and its environment.
5. **Interface CP Detection**: Assigns each CP to the interface between two fragments using a "nearest neighbor" geometric algorithm, without the need to set an arbitrary distance threshold.
6. **VMD Script Generation**: Produces a `.tcl` script per complex with all interface CPs pre-selected and highlighted. Each CP type is color-coded: **Silver** for Nuclear (3,−3), **Lime** for Bond (3,−1), **Purple** for Ring (3,+1), and **Yellow** for Cage (3,+3).
7. **Property Calculation**: Extracts 25 QTAIM topological properties at each interface CP.
8. **CSV Export**: Consolidates all per-CP property files into one summary `.csv` per complex.

---

## Calculated Properties

The following topological properties are extracted at each Critical Point and saved in `{name}_cp_properties.csv`:

- **Topological (QTAIM):**
  - Electron Density (ρ)
  - Laplacian of Electron Density (∇²ρ)
  - Lagrangian (*G*) and Hamiltonian (*K*) Kinetic Energy Densities
  - Potential Energy Density (*V*)
  - Energy Density (*H*)
  - Electron Localization Function (ELF) & Local Orbital Locator (LOL)
  - Average Local Ionization Energy (ALIE)
  - Electrostatic Potential (ESP)

- **NCI & Related:**
  - Reduced Density Gradient (RDG)
  - Sign(λ₂)·ρ (NCI coloring variable)
  - Interaction Region Indicator (IRI)
  - Delta-g (promolecular and Hirshfeld)

- **Additional:**
  - Alpha and Beta Electron Densities
  - Spin Density
  - Local Information Entropy
  - van der Waals Potential (C probe)
  - Wavefunction value for orbital 1

---

## Output Files

| File | Description |
|---|---|
| `output/{name}.pdb` | Molecular structure |
| `output/{name}_CPs.pdb` | All critical points |
| `output/{name}_paths.pdb` | Bond paths |
| `output/{name}.csv` | Atom coordinates table |
| `output/{name}_CPs.csv` | CP coordinates and types |
| `output/{name}_paths.csv` | Path point coordinates |
| `output/{name}_CPs.txt` | Raw CP data |
| `output/{name}.tcl` | VMD visualization script |
| `output/cp_properties/{name}_cp_{id}.txt` | Raw output per CP |
| `output/cp_properties/{name}_cp_properties.csv` | Summary table of all interface CP properties |

---

## Installation

### Option 1: Conda Environment (Recommended)

```bash
git clone https://github.com/1JELC1/Interaction-Map-Explorer.git
cd Interaction-Map-Explorer
conda env create -f environment.yml
conda activate ime-env
```

### Option 2: Manual Installation

- Python 3.9+
- Install dependencies:

```bash
pip install -r requirements.txt
```

> **Note:** `vedo` can occasionally have issues with certain PyQt backends. If you encounter display errors, the Conda environment is more reliable.

### External Requirement

- **Multiwfn**: This software cannot be installed automatically. Please download it from the [official website](http://sobereva.com/multiwfn) and place the Multiwfn (executable) in the same folder as the scripts or add it to your system `PATH`.

- **FragmentFinder**: Place `FragmentFinder.py` in the same folder as the scripts. Available at [1JELC1/FragmentFinder](https://github.com/1JELC1/FragmentFinder).

---

## Usage

### 1. Prepare your files

Create a folder (default: `wf/`) containing your wavefunction files:
```
wf/
├── complex_01.fchk
├── complex_02.fchk
└── ...
```

For **Batch Mode**, also prepare a single `.xyz` file of the reference ligand:
```
ligand.xyz
```

### 2. Run the script

```bash
python IME.py
```

### 3. Interactive Setup

The script will guide you through:

- **Directory paths**: Enter the wavefunction folder path and output folder name.
- **Mode selection**:
  - `[1] Interactive Mode` — Opens a 3D viewer for each file where you can manually define multiple fragments (minimum 2, up to 9 using the numeric keys `1`–`9`). The tool will search for CPs at every interface between the defined fragments.
  - `[2] Batch Mode` — Provide the reference ligand `.xyz` file once; the tool finds it automatically in all complexes and assigns it as **Fragment 1**. The remaining atoms are assigned as **Fragment 2** (environment). The tool then searches for interactions between the ligand and its environment.

- *(Batch Mode only)* **CP types**: Select which CP types to extract (Nuclear, Bond, Ring, Cage).

### 4. Interactive Mode controls (FragmentFinder)

In the interactive graphical interface, the user defines the fragments between which to search for interactions. The system is controlled via mouse clicks on atoms and keyboard shortcuts:

**Fragment Editing Mode:**

| Key | Action |
|---|---|
| Left-click | Select/deselect an atom for the current fragment |
| `1`–`9` | Switch to Fragment 1 through 9 |
| `n` | Add neighboring atoms to the current fragment |
| `m` | Clear all atoms from the current fragment |
| `k` | Assign all remaining unassigned atoms to the current fragment |
| `e` | Toggle atom labels |
| `f` | Run analysis: detect CPs and paths between the defined fragments |
| `q` | Finish and proceed to property calculation |
| `h` | Show help |

**CP Visualization Mode** (after pressing `f`):

| Key | Action |
|---|---|
| `z` | Toggle Nuclear CPs (3,−3) visibility |
| `x` | Toggle Bond CPs (3,−1) visibility |
| `c` | Toggle Ring CPs (3,+1) visibility |
| `v` | Toggle Cage CPs (3,+3) visibility |
| `d` | Enter/exit delete mode (click a CP to remove it) |
| `r` | Restore all previously deleted CPs |
| `f` | Return to fragment editing mode |
| `q` | Finish selection and proceed to property calculation |

### 5. Output

After processing, the `output/` folder will contain all PDB, CSV, TCL, and property files ready for VMD visualization and statistical analysis.

---

## Case Study: Protein-Ligand Interface

Here is a demonstration of how to define fragments in **Interactive Mode** using a ligand inside a protein active site.

**Example 1: Ligand vs. Entire Active Site**
The ligand is selected as Fragment 1, and the rest of the active site atoms are assigned to Fragment 2. IME finds all the interaction CPs between the ligand and the pocket.

https://github.com/user-attachments/assets/55121fb1-85db-492b-b500-23a7e93fe7f5

**Example 2: Triad Interaction (Ligand and Two Specific Residues)**
The ligand is Fragment 1, and two specific surrounding residues are individually selected as Fragment 2 and Fragment 3. This isolates and characterizes only the interactions forming this specific triad.

https://github.com/user-attachments/assets/d8dec599-0630-4281-b449-4dface673289

---

## VMD Visualization

Open VMD and load the generated `.tcl` script:

```tcl
source output/complex_01.tcl
```

The script automatically loads the molecule, all CPs (hidden by default), and highlights the **interface CPs** with color-coded representations:

| CP Type | VMD Color | Description |
|---|---|---|
| (3,−3) Nuclear | Silver | Atomic critical points |
| (3,−1) Bond | Lime | Bond critical points |
| (3,+1) Ring | Purple | Ring critical points |
| (3,+3) Cage | Yellow | Cage critical points |

---

## Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](https://github.com/1JELC1/Interaction-Map-Explorer/issues) if you want to report a bug or request a feature.

---

## Citation

If you use this software in your research, please cite it using our Zenodo DOI:

**DOI**: [10.5281/zenodo.19100311](https://doi.org/10.5281/zenodo.19100311)

You can also use the citation metadata provided in the [CITATION.cff](CITATION.cff) file.

---

## License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](LICENSE) file for details.
