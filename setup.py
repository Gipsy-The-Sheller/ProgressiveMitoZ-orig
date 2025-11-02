#!/usr/bin/env python3
from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()


setup(
    name='pmitoz',
    version = '3.6',
    author = 'Zhi-Jie Xu',
    author_email = 'zjxmolls@outlook.com',
    long_description = long_description,
    license = 'GPLv3',
    packages=find_packages(),
    entry_points = {'console_scripts': ['pmitoz = pmitoz.MitoZ:main', 'mitoz-tools = pmitoz.tools:main']},
    url='https://github.com/linzhi2013/MitoZ',
    python_requires='>=3',
    install_requires = [
	    # No dependencies listed here since we need to rely on conda anyway
    ],
    include_package_data=True, # https://setuptools.pypa.io/en/latest/userguide/datafiles.html. This tells setuptools to install any data files it finds in your packages. The data files must be specified via the MANIFEST.in file.
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPLv3 License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)


