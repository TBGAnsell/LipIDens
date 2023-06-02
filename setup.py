"""
LipIDens

A protocol for simulation assissted interpretation of lipid-like densities.

"""

from setuptools import setup, find_packages
import re

with open("README.md", "r") as fle:
    ld=fle.read()

try:
    with open("lipidens/_version.py", "rt") as v_fle:
        v_str=v_fle.read()
        v_s=r"^__version__ = ['\"]([^'\"]*)['\"]"
        v=re.search(v_s, v_str, re.M)
        curr_ver=v.group(1)
except Exception as e:
    raise RuntimeError("Unable to locate version string:", e)

import os

def get_files(direc):
    fle_list=[]
    for (path, _, fle_name) in os.walk(direc):
        for f_n in fle_name:
            fle_list.append(os.path.join('..', path, f_n))
    return fle_list


setup(
  name = 'lipidens',
  version = curr_ver,
  license="MIT",
  author = 'T. Bertie Ansell, Wanling Song',
  author_email = 'bertie.ansell@bioch.ox.ac.uk',
  url = 'https://github.com/TBGAnsell/LipIDens',
  description="A protocol for simulation assisted interpretation of lipid-like densities.",
  long_description=ld,
  long_description_content_type="text/markdown",
  packages=find_packages(),
  include_package_data=True,
  package_data={"": get_files('lipidens')},
  keywords = ['simulation', 'lipid', 'density', 'binding site'],
  python_requires='>=3.9, <4',
  install_requires=['kneebow',
      'logomaker',
      'matplotlib>=3.3.4',
      'mdtraj',
      'networkx',
      'numpy',
      'p-tqdm',
      'pandas',
      'pylipid>=1.5.14',
      'python-louvain',
      'scikit-learn',
      'scipy',
      'seaborn',
      'statsmodels',
      'tqdm',
      'vermouth>=0.7.2'
      ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Medical Science Apps.',
    'Topic :: Scientific/Engineering :: Mathematics',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.9',
  ],
)
