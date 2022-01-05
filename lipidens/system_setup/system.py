#!/usr/bin/env python3

import os

"""
Basic function for establishing path of simulation data.
"""

def setup(protocol_path, save_dir=None):
    """
    Establish save directory for simulations and data.
    """
    if save_dir is None:
        path=os.getcwd()
    else:
        if os.path.isabs(save_dir)==True:
            path=save_dir
        else:
            path=os.path.join(os.getcwd(), save_dir)

    print("Location of data:", path)
    return path
