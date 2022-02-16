from setuptools import setup, find_packages

setup(
    author = "Tobias von der Haar",
    description = "A Python package for gene optimisation",
    name = "CodOpY",
    version = "0.2.0",
    packages = find_packages(include = ["CodOpY","CodOpY.*"]),
    include_package_data = True,
    install_requires = ['pandas','numpy','matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)