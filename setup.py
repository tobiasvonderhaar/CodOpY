from setuptools import setup, find_packages


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
print(find_packages())
setup(
    name = 'CodOpY',
    version = '0.1.1',
    author = 'Tobias von der Haar',
    author_email = 'T.von-der-Haar@kent.ac.uk',
    url="https://github.com/tobiasvonderhaar/CodOpY",
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
    py_modules = ['CodOpY.optimise','CodOpY.misc','CodOpY.plot'],
    zip_safe=False,
    package_dir = {'' : 'src'},
    packages=find_packages(include=['CodOpY','CodOpY.*']),
    include_package_data=True,
    package_data={'':['src/CodOpY/Data/*.csv']},
    python_requires = ">=3.6",
    install_requires(['pandas','numpy','matplotlib'])
)
