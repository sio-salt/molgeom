# molgeom ⌬  
A Simple Python library for molecular geometry manipulation.


## ✨ Key Features
- Read and write molecular geometries:
  - Formats supported: XYZ, Gaussian input, GAMESS input, CIF, POSCAR.
- Calculate molecular properties:
  - Nuclear repulsion energy.
  - Center of mass.
  - Bond lengths.
- Perform geometry manipulations:
  - Translation and rotation (including matrix and Rodrigues' rotation formula).
  - Cell replication using lattice vectors and symmetry operations (e.g., `-y, x + 1/2, -z + 1/2`).
  - Merge molecule geometries.
- Additional utilities:
  - Filter atoms by element.
  - Cluster fragments.
  - Identify bond cycles.


## 🚀 Getting Started
You can try out example of `molgeom` in Jupyter Notebook:  
👉 [**Try Example Code in Jupyter Notebook**](https://mybinder.org/v2/gh/sio-salt/molgeom-examples/main?labpath=notebooks%2Ftutorial1.ipynb)

For a version with NumPy, see the [**numpyfy branch**](https://github.com/sio-salt/molgeom/tree/numpyfy).  

---

## Installation
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

## Future Plans
- Support Z-matrix format.
- Add more symmetry utilities.
- Handle multi-molecule systems.
- Additional features and improvements.

## Links
- [GitHub Repository](https://github.com/sio-salt/molgeom/tree/main)
- [Numpyfy Branch](https://github.com/sio-salt/molgeom/tree/numpyfy)
- [Example Repository](https://github.com/sio-salt/molgeom-examples/tree/main)
- [Binder Link](https://mybinder.org/v2/gh/sio-salt/molgeom-examples/main?labpath=notebooks%2Ftutorial1.ipynb)
