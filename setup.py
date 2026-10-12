#!/usr/bin/env python3
"""
Setup script for Echidna
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read version from __init__.py
import re
_version_match = re.search(
    r'^__version__\s*=\s*["\']([^"\']+)["\']',
    open("echidna/__init__.py").read(),
)
version = {"__version__": _version_match.group(1) if _version_match else "0.0.0"}

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = ""
if readme_file.exists():
    long_description = readme_file.read_text()

setup(
    name="echidna",
    version=version.get("__version__", "1.0.0"),
    description="Echidna - metagenomic assembly, MAG integrity and chimera screening",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Luigi Gallucci",
    author_email="luigigallucci8@gmail.com",
    url="https://github.com/l-gallucci/echidna",
    license="MIT",

    packages=find_packages(),

    python_requires=">=3.8",
    
    install_requires=[
        # No external dependencies required for core functionality
    ],
    
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.10',
            'black>=21.0',
            'flake8>=3.8',
        ],
    },
    
    entry_points={
        'console_scripts': [
            'echidna=echidna.__main__:main',
        ],
    },
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],

    keywords="bioinformatics metagenomics assembly quality-control",

    project_urls={
        "Bug Reports": "https://github.com/l-gallucci/echidna/issues",
        "Source": "https://github.com/l-gallucci/echidna",
    },
)
