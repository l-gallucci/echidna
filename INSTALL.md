# Installation Guide

## Quick Install

```bash
# Clone repository
git clone https://github.com/l-gallucci/echidna.git
cd echidna

# Install in development mode
pip install -e .

# Or install normally
pip install .
```

## Verify Installation

```bash
echidna --version
echidna --help
```

## Directory Structure

After cloning, your directory should look like:

```
echidna/
├── echidna/
│   ├── __init__.py
│   ├── __main__.py
│   ├── commands/
│   │   ├── __init__.py
│   │   ├── filter_misassemblies.py
│   │   ├── remove_hq_contigs.py
│   │   ├── filter_checkm2.py
│   │   └── run_gtdbtk.py
│   └── utils/
│       ├── __init__.py
│       ├── fasta.py
│       ├── logging.py
│       └── validation.py
├── setup.py
├── README.md
├── LICENSE
└── INSTALL.md
```

## Dependencies

### Core (required)
- Python >= 3.7

### External Tools (optional, depending on commands used)

#### For `filter-misassemblies`:
```bash
# Install anvi'o
conda create -n anvio-8 -c conda-forge -c bioconda anvio=8
```

#### For `filter-checkm2`:
```bash
# Install CheckM2
conda create -n checkm2 -c conda-forge -c bioconda checkm2
```

#### For `run-gtdbtk`:
```bash
# Install GTDB-Tk
conda create -n gtdbtk-2 -c conda-forge -c bioconda gtdbtk=2.3.2

# Download GTDB-Tk database
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz
export GTDBTK_DATA_PATH=/path/to/gtdbtk_data
```

## Development Installation

If you want to contribute to Echidna:

```bash
# Clone repository
git clone https://github.com/l-gallucci/echidna.git
cd echidna

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode with development dependencies
pip install -e ".[dev]"

# Run tests (when available)
pytest
```

## Troubleshooting

### Command not found after installation

If `echidna` is not found after installation:

```bash
# Check if it's in your PATH
which echidna

# If not, find where pip installed it
pip show echidna

# Add to PATH or use full path
export PATH="$PATH:/path/to/echidna/bin"
```

### Permission denied

If you get permission errors:

```bash
# Install for user only
pip install --user .

# Or use sudo (not recommended)
sudo pip install .
```

### Python version issues

Check your Python version:

```bash
python --version  # Should be >= 3.7

# If you have multiple Python versions
python3.9 -m pip install .
```

## Uninstallation

```bash
pip uninstall echidna
```

## Platform-Specific Notes

### Linux
No special requirements.

### macOS
Works out of the box. If using conda, prefer:
```bash
conda install -c conda-forge python=3.9
```

### Windows
Echidna should work on Windows, but external tools (anvi'o, CheckM2, GTDB-Tk) are best run in WSL2 or Docker.

## Docker (Coming Soon)

```bash
# Build Docker image
docker build -t echidna:latest .

# Run
docker run -v $(pwd):/data echidna filter-misassemblies --help
```
