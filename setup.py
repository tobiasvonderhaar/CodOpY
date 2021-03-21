from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name = 'CodOpY',
    version = '0.0.1',
    author = 'Tobias von der Haar'
    author_email = 'T.von-der-Haar@kent.ac.uk'
    url="https://github.com/tobiasvonderhaar",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],

    description= 'Python Codon Optimisation tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    py_modules = ['CodOpY.optimise','CodOpY.misc'],
    package_dir = {'' : 'src'},
    python_requires = ">=3.6"
)
