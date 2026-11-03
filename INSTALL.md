# Installation

The Python core requires Python >=3.8 and no third-party runtime libraries.
External bioinformatics tools and their databases are separate dependencies;
choose an interpreter supported by each tool's environment.

```bash
git clone https://github.com/l-gallucci/echidna.git
cd echidna
python -m pip install -e .
echidna --version
echidna --help
```

See [the installation guide](docs/installation.md) for tool environments and
reference databases, and [the wiki](wiki/Home.md) for workflow contracts.

The distribution is `echidnaxomics`; the command and Python import are `echidna`.
Until publication, install from this checkout. Do not install the unrelated
`echidna` distribution from a registry.

Wheels include Snakemake templates. Export a writable copy after installation:

```bash
echidna setup-workflow --export-workflow ./my-workflow
```

Use a new destination; export never overwrites an existing directory. With a
checkout, the existing `--repo /absolute/path/to/echidna` mode still works and
points the tool environments at that editable checkout. Configure exported
`config.yaml`, then run Snakemake with its absolute `-s` and `--configfile` paths.
Input/output paths remain relative to the directory where Snakemake is launched.
Pass `fetch-databases --config /absolute/path/to/my-workflow/config.yaml` when
registering databases for an exported workflow.
See [release preparation](docs/releasing.md) for author metadata and Bioconda.

For development, use an environment with pytest:

```bash
python -m pip install -e ".[dev]"
PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider
```

Python core portability does not guarantee availability of every external tool
on macOS or Windows. Follow upstream requirements; use Linux/WSL or suitable
containers where needed. The project currently ships no Dockerfile.
