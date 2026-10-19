# LipIDens - simulation assisted interpretation of lipid densities 

![logo](lipidens_final_Asset_2.png)

**LipIDens** is an open-source pipeline for **simulation-assisted interpretation of lipid and lipid-like densities** observed in cryo-electron microscopy (cryo-EM) structures of membrane proteins.

LipIDens integrates **coarse-grained and atomistic molecular dynamics (MD) simulations** with **PyLipID-based analysis** to identify, characterise, and rank lipid binding sites, enabling experimental densities to be connected to specific lipid chemistries.

---

## Key Capabilities

- Differentiating between sterol and phospholipid binding sites where structural densities are ambiguous. 
- Assess preferential binding of different lipid types to a site where lipid-like density is observed.
- Evaluate whether lipid tail densities are likely to belong to the same or adjacent binding sites. 
- Quantify the kinetics of different lipids binding to the same site.
- Obtain a more complete picture of lipid interactions profiles around the protein in a membrane environment. 
- Assist interpretation of cryo-EM densities during model building and map refinement cycles. 
- Assess which putative lipid binding sites may prevail when not in detergent conditions. 
---

## Scope of the Pipeline

LipIDens provides **workflow automation and analysis tools**. The user is responsible for:

- Running MD simulations (e.g. with **GROMACS**)
- Providing trajectories and topology files

LipIDens then processes these data to guide structural interpretation.

---

## Installation

Python **≥ 3.11** is required.

```bash
mamba create -n lipidens python=3.12 -y
mamba activate lipidens
mamba install -c conda-forge gromacs pdb2pqr -y
pip install git+https://github.com/pstansfeld/LipIDens.git
```

This installs LipIDens with its Python dependencies, including [PyLipID](https://github.com/pstansfeld/PyLipID) (≥ 1.6)
and [vermouth](https://github.com/marrink-lab/vermouth-martinize) (`martinize2`).

For the atomistic stages (backmapping lipid poses and setting up atomistic simulations) also install
[cg2at-lite](https://github.com/pstansfeld/cg2at-lite) (`pip install cg2at-lite`); the original
[CG2AT](https://github.com/owenvickery/cg2at) is used instead if cg2at-lite is not installed. Comparing poses with
densities (stage 5) writes a [PyMOL](https://pymol.org) script.

Optional:
```bash
mamba install jupyter -y
```

---

## Running LipIDens

Run the protocol interactively, stage by stage:

```bash
lipidens
```

(or `python lipidens_master_run.py` from a clone of this repository), or follow the
[LipIDens workflow notebook](https://github.com/pstansfeld/LipIDens/blob/main/LipIDens.ipynb).

If you write your own scripts that call LipIDens or PyLipID, put the code under `if __name__ == "__main__":`:
PyLipID calculates in worker processes that re-import the script.

---

## Typical Workflow

1. Prepare a membrane protein system and run coarse‑grained MD simulations
2. Analyse lipid interactions using PyLipID via LipIDens
3. Screen and rank lipid binding sites
4. Refine selected lipid poses using atomistic simulations
5. Interpret lipid‑like cryo‑EM densities in a bilayer context

A worked example of the full LipIDens protocol is provided in the **[LipIDens workflow notebook](https://github.com/pstansfeld/LipIDens/blob/main/LipIDens.ipynb)**.

---

## Dependencies

- Python ≥ 3.11
- PyLipID ≥ 1.6, mdtraj, numpy, pandas, matplotlib, seaborn, vermouth (installed by pip)
- GROMACS and pdb2pqr (user-installed; simulations are run by the user)
- cg2at-lite or CG2AT (atomistic stages)
- MARTINI force field (downloaded by LipIDens)

---

## Citation

**Please cite** the following if you use LipIDens in your research:

Ansell, T.B., Song, W., Coupland, C.E., Carrique, L., Corey, R.A., Duncan, A.L., Cassidy, C.K., Geurts, M.M.G., Rasmussen, T., Ward, A.B., Siebold, C., Stansfeld, P., Sansom, M.S.P. (2022). **LipIDens: Simulation assisted interpretation of lipid densities in cryo-EM structures of membrane proteins.** *Nature Communications* **14**, 7774 (2023).  [doi: 10.1038/s41467-023-43392-y](https://www.nature.com/articles/s41467-023-43392-y)

**Accompanying step by step protocol** citation:

Ansell, T.B., Song, W., Coupland, C.E., Carrique, L., Corey, R.A., Duncan, A.L., Cassidy, C.K., Geurts, M.M.G., Rasmussen, T., Ward, A.B., Siebold, C., Stansfeld, P., Sansom, M.S.P. (2023). **Implementation of the LipIDens pipeline: assisted interpretation of lipid densities in membrane protein structures using simulations.** protocols exchange, [doi: 10.21203/rs.3.pex-2408/v1](https://doi.org/10.21203/rs.3.pex-2408/v1)

The LipIDens logo was designed by Jessica Ansell. 

---

## License

LipIDens is released under the MIT License.
