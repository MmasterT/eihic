#!/usr/bin/env python

from setuptools import setup
import glob

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]

setup(
    name="eihic",
    version="0.1.0",
    description=" EI Hi-C QC pipeline",
    author="Mariano Olivera Fedi",
    author_email="Mariano.olivera-fedi@aerlham.ac.uk, mariano.olivera.fedi@hotmail.com",
    url="https://github.com/EI-CoreBioinformatics/eihic",
    license="CC BY-NC 4.0",
    zip_safe=False,
    keywords="Hi-C Arima Dovetail omni-c",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.9",
    ],
    entry_points={"console_scripts": ["eihic = eihic.__main__:main"]},
    install_requires=requirements,
    packages=["eihic", "eihic.etc", "eihic.scripts", "eihic.workflow"],
    scripts=[script for script in glob.glob("eihic/scripts/*")],
    package_data={
        "eihic.workflow": ["omni-c.smk",
            "arima.smk"],
        "eihic.etc": ["hpc_config.json", 
            "run_config.yaml", 
            "samples.csv"],
    },
    include_package_data=True,
)