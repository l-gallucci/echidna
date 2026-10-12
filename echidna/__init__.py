"""
Echidna - metagenomic assembly, MAG integrity and chimera screening
A toolkit for metagenomic assembly quality control and filtering
"""

__version__ = "2.0.0"
__author__ = "Luigi Gallucci"
__license__ = "MIT"

from pathlib import Path

# Package root directory
PACKAGE_ROOT = Path(__file__).parent

__all__ = ['__version__', '__author__', '__license__', 'PACKAGE_ROOT']
