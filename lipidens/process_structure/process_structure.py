#!/usr/bin/env python3

import subprocess
import sys
import os

def process_structure(protocol_path, protein_AT_full, protein_rotate, path):
    """
    Call the bash script process.sh to use GROMACS to orient the protein along principle axis. Returns orientated structure.

    Params:
    -------
    protocol_path: str
        path to LipIDens
    protein_AT_full: str
        input atomistic protein coordinate file
    protein_rotate: str
        degrees to rotate protein in 'x y z'
    path: str
        path
    
    Returns:
    --------
    protein_AT_full: str
        oriented protein coordinate file

    """
    if os.path.isfile(f"{protein_AT_full}"):
            protein_AT_full=subprocess.check_output([f"{protocol_path}/process_structure/process.sh", protein_AT_full, protein_rotate]).decode(sys.stdout.encoding).strip()
    else:
        print(f"{protein_AT_full} not found - check protein .pdb file exists and in correct location")
        quit()
    return protein_AT_full
