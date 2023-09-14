# setup.py

from setuptools import setup, find_packages

setup(
    name='Bio.Informatica',
    version='1.0',
    packages=find_packages(),
    install_requires=[
    'biopython',
    'neatbio',  # Assuming this is a valid package name
    'numpy',
    'pandas',
    'rdkit',
    'altair',  # Optional, include if needed
    'mols2grid',  # Assuming this is a valid package name
],
)
