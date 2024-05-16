from setuptools import setup, find_packages

setup(
    name="rna3db",
    version=1.1,
    description="A dataset for training and benchmarking deep learning models for RNA structure prediction",
    author="Marcell Szikszai",
    packages=find_packages(exclude=["tests", "scripts", "data"]),
    install_requires=["biopython", "tqdm", "black", "pre-commit"],
)
