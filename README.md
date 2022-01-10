# LipIDens - simulation assisted interpretation of lipid densities 

LipIDens is a pipeline for simulation assisted interpretation of lipid or lipid-like densities in e.g. cryogenic electron microscopy (cryo-EM) structures of membrane proteins. The protocol can be used to establish CG simulations, analyse lipid binding sites using PyLipID, screen, rank and process PyLipID outputs and establish atomistic simulations to refine lipid binding poses.

## Applications:
- Differentiating between sterol and phospholipid binding sites where structural densities are ambiguous. 
- Assess preferential binding of different lipid types to a site where lipid-like density is observed.
- Evaluate whether lipid tail densities are likely to belong to the same or adjacent binding sites. 
- Quantify the kinetics of different lipids binding to the same site.
- Obtain a more complete picture of lipid interactions profiles around the protein in a membrane environment. 
- Assist interpretation of cryo-EM densities during model building and map refinement cycles. 
- Assess which putative lipid binding sites may prevail when not in detergent conditions. 

## Installation:

LipIDens requires a python3 environment (>=3.9 recommended). 
```bash
pip install lipidens
```
```bash
git clone XXX
```
```bash
conda create -n LipIDens python=3.9
pip install -r requirements.txt
```
## Usage

Detailed steps for the usage and implementation of LidIDens are provided within the accompanying protocols manuscript (see citation below). 

In addition, a jupyter notebook tutorial is provided to assist implementation. 

LipIDens can also be run using the standalone master_run.py script.  


## Citation:

Citation to follow. 

