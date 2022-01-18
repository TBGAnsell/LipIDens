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

def run_type():
    steps_d={"1a": "Setup coarse-grained simulations",
             "1b": "Process coarse-grained trajectories",
             "2": "Test PyLipID cut-off values",
             "3": "Run PyLipID analysis",
             "4": "Screen PyLipID data",
             "5": "Rank site lipids",
             "6": "Setup atomistic simulations"}
    print("\nInitialising LipIDens protocol\nProtocol steps:")
    print(*[f"\n{key}: {val}" for key, val in steps_d.items()])
    try:
        input_step=str(input("\nSelect protocol stage number:"))
        if input_step in steps_d:
            print("\nSelected:"+steps_d[input_step])
        else:
            print("Invalid: enter number of required protocol stage")
            exit()
    except Exception as e:
        print(e)
        exit()
    return input_step

