# molgeom ⌬
A Simple Python library for molecular geometry manipulation.

> **Note**: `molgeom` is under active development. The implemented features and method names may change in future updates. Please check the repository regularly for the latest information.

![molgeom_view_mols_example_1](https://github.com/user-attachments/assets/c6e7775c-6e07-4c99-8760-e4f7b7cc2679)


## ✨ Key Features
- Read and write molecular geometries:
  - Formats supported: XYZ, Gaussian input, GAMESS input, MOL, SDF, CIF, VASP POSCAR.
- Calculate molecular properties:
  - Nuclear repulsion energy.
  - Center of mass.
  - Bond lengths.
- Perform geometry manipulations:
  - Translation and rotation (including matrix and Rodrigues' rotation formula).
  - Cell replication using lattice vectors and symmetry operations (e.g., `-y, x + 1/2, -z + 1/2`).
  - Merge molecule geometries.
- Visualize molecular geometries:
  - Display molecules in your Jupyter Notebook or browser using 3Dmol.js.
- Command-line interface:
  - `molgeom` command is available after installation.
  - See Command Line Usage section below for details
- Additional utilities:
  - Filter atoms by element.
  - Cluster fragments.
  - Identify bond cycles.



## 🚀 Getting Started

You can try out example of `molgeom` in Jupyter Notebook:
👉 [**Try Example Code in Jupyter Notebook**](https://mybinder.org/v2/gh/sio-salt/molgeom-examples/main?urlpath=lab/tree/notebooks/tutorial1.ipynb)

For a version with NumPy, see the [**numpyfy branch**](https://github.com/sio-salt/molgeom/tree/numpyfy).

---

## 🔽 Installation
Install `molgeom` using pip:
```bash
pip install -U git+https://github.com/sio-salt/molgeom@main
```
Alternatively, clone the repository for development:

```bash
git clone https://github.com/sio-salt/molgeom.git
cd molgeom
pip install -e .
```
The minimum Python version required is 3.9 for the main branch.

## 📟 Command Line Usage
After installation, the `molgeom` command becomes available. Enable shell completion by copying and pasting the following command into your terminal:

Bash:
```bash
echo 'eval "$(_MOLGEOM_COMPLETE=bash_source molgeom)"' >> ~/.bashrc && source ~/.bashrc
```

Zsh:
```zsh
echo 'eval "$(_MOLGEOM_COMPLETE=zsh_source molgeom)"' >> ~/.zshrc && source ~/.zshrc
```

Fish:
```fish
_MOLGEOM_COMPLETE=fish_source molgeom | source
```

Available commands:
```bash
# Single file commands
molgeom modify <file> [-op <op>]   # Transform structure (translate/reflect/rotate)
molgeom split <file>               # Split into molecular clusters
molgeom poscar2xyz <file> <ranges> # Convert POSCAR with cell replication
                                   # ( e.g., molgeom poscar2xyz POSCAR -1 2 -1 2 -1 2 )

# Multiple file commands (space-separated)
molgeom center <file1> <file2> ...         # Print center of mass for multiple files
molgeom nuclrep <file1> <file2> ...        # Calculate nuclear repulsion energy for multiple files
molgeom bonds <file1> <file2> ... [--tol]  # List bonds for multiple files
molgeom view <file1> <file2> ...           # View multiple structures in browser

# Get help for any command
molgeom --help
molgeom <command> --help
```

## 🛠️ Future Plans
- Support Z-matrix format.
- Add more symmetry utilities.
- Handle multi-molecule systems.
- Additional features and improvements.


## 🔗 Links
- [main branch](https://github.com/sio-salt/molgeom/tree/main)
- [numpyfy branch](https://github.com/sio-salt/molgeom/tree/numpyfy)
- [molgeom examples repo](https://github.com/sio-salt/molgeom-examples/tree/main)
- [molgeom examples Binder](https://mybinder.org/v2/gh/sio-salt/molgeom-examples/main?urlpath=lab/tree/notebooks/tutorial1.ipynb)
