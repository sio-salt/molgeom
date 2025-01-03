# molgeom ⌬  
A simple and powerful Python library for molecular geometry manipulation.  

## 🌟 Key Features
- **Read and write molecular geometries** from popular file formats:
  - XYZ, Gaussian input, GAMESS input, CIF, POSCAR.
- **Manipulate geometries**:
  - Translation and rotation (supports matrix and Rodrigues' rotation).
  - Cell replication using lattice vectors or symmetry operations (e.g., `-y, x + 1/2, -z + 1/2`).
  - Merging molecule geometries.
- **Calculate molecular properties** with ease:
  - Nuclear repulsion energy.
  - Center of mass.
  - Bond detection.
- **Additional utilities**:
  - Filter atoms by element.
  - Cluster fragments.
  - Identify bond cycles.


## 🚀 Get Started
Want to try `molgeom` now? Explore the tutorial directly in your browser:  
👉 [**Try Example Code in Jupyter Notebook**](https://mybinder.org/v2/gh/sio-salt/molgeom-examples/main?labpath=notebooks%2Ftutorial1.ipynb)


For advanced performance, check out the [**numpyfy branch**](https://github.com/sio-salt/molgeom/tree/numpyfy).

---

## 📥 Installation
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
- Add more symmetry operations.
- Handle multi-molecule systems.
- ...and much more.

## 🔗 Links
- [GitHub Repository](https://github.com/sio-salt/molgeom)
- [Numpyfy Branch](https://github.com/sio-salt/molgeom/tree/numpyfy)
- [Example Code](https://github.com/sio-salt/molgeom-examples/tree/main)
