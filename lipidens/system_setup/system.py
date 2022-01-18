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
    steps_d={1: "Setup coarse-grained simulations",
             2: "Process coarse-grained trajectories",
             3: "Analyse simulations",
             4: "Setup atomistic simulations"}
    print("\nInitialising LipIDens protocol\nProtocol steps:")
    print(*[f"\n{key}: {val}\n" for key, val in steps_d.items()])
    try:
        input_step=int(input("Select protocol stage number:"))
        if input_step in steps_d:
            print("\nSelected:"+steps_d[input_step])
        else:
            print("Invalid: enter number of required protocol stage")
            exit()
    except Exception as e:
        print(e)
        exit()
    return input_step

