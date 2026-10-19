## Installation

LipIDens requires **Python ≥ 3.11**.

Clone and install with pip:

```bash
git clone https://github.com/pstansfeld/LipIDens.git
cd LipIDens
pip install .
```

For development (editable install):

```bash
pip install -e .
```

Either installs the `lipidens` command, which runs the protocol interactively, together with the workflow scripts and
parameter files it uses. GROMACS, pdb2pqr and cg2at-lite (or CG2AT) are installed separately: see README.md.
**MD simulations (e.g. GROMACS) must be run separately by the user.**
