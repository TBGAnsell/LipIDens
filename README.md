# LipIDens - simulation assisted interpretation of lipid densities 

![logo](lipidens_final_Asset_2.png)

LipIDens is an open-source pipeline for simulation-assisted interpretation of lipid and lipid-like densities observed in cryo-electron microscopy (cryo-EM) structures of membrane proteins.

LipIDens integrates coarse-grained and atomistic molecular dynamics (MD) simulations with PyLipID-based analysis to identify, characterise, and rank lipid binding sites, enabling experimental densities to be connected to specific lipid chemistries. 

## Key Capabilities:
- Differentiating between sterol and phospholipid binding sites where structural densities are ambiguous. 
- Assess preferential binding of different lipid types to a site where lipid-like density is observed.
- Evaluate whether lipid tail densities are likely to belong to the same or adjacent binding sites. 
- Quantify the kinetics of different lipids binding to the same site.
- Obtain a more complete picture of lipid interactions profiles around the protein in a membrane environment. 
- Assist interpretation of cryo-EM densities during model building and map refinement cycles. 
- Assess which putative lipid binding sites may prevail when not in detergent conditions. 

## Scope of the Pipeline:
LipIDens provides workflow automation and analysis tools. The user is responsible for:

Running MD simulations (e.g. with GROMACS)
Providing trajectories and topology files
LipIDens then processes these data to guide structural interpretation.

## Installation:

For detailed instructions on LipIDens Installation please refer to the latest version, courtesy of the Stansfeld Lab, University of Warwick.

[Latest LipIDens Version and Installation](https://github.com/pstansfeld/LipIDens/), updated April 2026

LipIDens requires a python3 environment (>=3.9 recommended). 

## Citation:

**Please cite** the following if you use LipIDens in your research:

Ansell, T.B., Song, W., Coupland, C.E., Carrique, L., Corey, R.A., Duncan, A.L., Cassidy, C.K., Geurts, M.M.G., Rasmussen, T., Ward, A.B., Siebold, C., Stansfeld, P., Sansom, M.S.P. (2022). **LipIDens: Simulation assisted interpretation of lipid densities in cryo-EM structures of membrane proteins.** Nature Communications 14, 7774, [doi: 10.1038/s41467-023-43392-y](https://www.nature.com/articles/s41467-023-43392-y)

**Accompanying step by step protocol** citation:

Ansell, T.B., Song, W., Coupland, C.E., Carrique, L., Corey, R.A., Duncan, A.L., Cassidy, C.K., Geurts, M.M.G., Rasmussen, T., Ward, A.B., Siebold, C., Stansfeld, P., Sansom, M.S.P. (2023). **Implementation of the LipIDens pipeline: assisted interpretation of lipid densities in membrane protein structures using simulations.** protocols exchange, [doi: 10.21203/rs.3.pex-2408/v1](https://doi.org/10.21203/rs.3.pex-2408/v1)

The LipIDens logo was designed by Jessica Ansell. 

## Legacy Installation (OLD):

See above for latest installation of LipIDens. 

LipIDens installation from the GitHub repository:
```bash
git clone https://github.com/TBGAnsell/LipIDens
cd LipIDens
python setup.py install
```

#### Recomended

Recommended to use dssp version 3 as this is the only current version compatible with martinize2. Check your version using `mkdssp --version`.

```bash
conda install -c salilab dssp
```


Conda (legacy installation):
```bash
conda create -n LipIDens python=3.9
pip install -r requirements.txt
```

