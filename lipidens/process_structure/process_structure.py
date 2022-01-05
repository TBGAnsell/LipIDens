#!/usr/bin/env python3

import subprocess
import sys

def process_structure(protocol_path, protein_AT_full, protein_rotate):
    """
    Call the bash script process.sh to use GROMACS to orient the protein along principle axis. Returns orientated structure.
    """
    protein_AT_full=subprocess.check_output(["{}/process_structure/process.sh".format(protocol_path), protein_AT_full, protein_rotate]).decode(sys.stdout.encoding).strip()
    return protein_AT_full
