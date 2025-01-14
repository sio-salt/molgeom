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
  - View molecules in 3D with 3Dmol.js in Jupyter Notebook or Browser.
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
