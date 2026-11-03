# Release preparation: echidnaxomics

This is a preparation guide. No upload, DOI, release, or Bioconda submission has
been made by the audit. The repository is still `l-gallucci/echidna`.

| Item | Name |
|---|---|
| Python/PyPI/Bioconda distribution | `echidnaxomics` |
| Shell command | `echidna` |
| Python import and `python -m` module | `echidna` |

`echiDNA` would normalize to `echidna`, a name already used by unrelated packages.
Do not install the unrelated `echidna` package from PyPI or Bioconda for this project.
The proposed name returned HTTP 404 from the PyPI JSON and Bioconda package APIs
on 2026-11-03; this does not reserve it or guarantee registry acceptance.

## Author metadata you should complete

Edit [`CITATION.cff`](../CITATION.cff): Luigi Gallucci is the sole author currently
listed. Add your real ORCID and affiliation at the commented locations, and any
coauthors as additional `authors` entries. Before a release, add `version` and
`date-released` matching its tag. The software version lives in
[`echidna/__init__.py`](../echidna/__init__.py). Current 0.0.4 is retained for local
validation; choose the next release version before uploading any artifacts.
Check the contact email in [`setup.py`](../setup.py). Add a real DOI to the CFF
only after obtaining it. A separate `.zenodo.json` is intentionally not supplied:
CFF is the single metadata source.

## Checks before releasing

The CI workflow runs Python tests, Snakemake integration fixtures, CFF validation,
source/wheel builds and installation outside the checkout. It does not publish.
Its remote result is only available after you push. A green package test does not
validate scientific performance or database compatibility.

```bash
python -m pip install build twine pytest cffconvert
python -m pytest -o addopts='' -q
cffconvert --validate
python -m build
python -m twine check dist/*
```

Use a fresh build directory after choosing a version; do not upload stale files
from another version. The core has no third-party runtime Python imports.
CheckM2, GTDB-Tk, GUNC, anvi'o and other subprocess tools have their own Python
and native dependencies (including pandas, numpy or scikit-learn depending on
the tool), installed in their separate environments. The base conda recipe does
not combine these incompatible tool stacks.

Before claiming the complete pipelines are validated, run both branches on
small, documented real datasets with the intended tools and databases, record
versions and checksums, and inspect clean and deliberately problematic controls.
The fixture DAG validates file handoffs, paths, optional RNA, tier changes and
reruns; it simulates expensive upstream results. It is not a real-database E2E
benchmark. The supplied Snakefile has a narrower scope than `mag-quality-wf`:
no MMseqs, SSU or final-report and no chimera-based exclusion.

## DOI: Zenodo software archive

When you decide to publish, connect the GitHub repository in Zenodo and enable
archiving **before** creating the intended GitHub release. Review the CFF metadata
first. A new GitHub release triggers an archive and DOI when that integration is
enabled. Zenodo supports both a version DOI and an all-versions DOI; cite the
version used for reproducibility. Review the archived record and add the actual
DOI to the repository metadata afterward.

An archive DOI identifies the software; it does not constitute peer review or
methodological validation. A methods article or protocol is a separate output.
No Zenodo integration or release has been activated here.

References: [Zenodo: CFF metadata](https://help.zenodo.org/docs/github/describe-software/citation-file/),
[archive a GitHub release](https://help.zenodo.org/docs/github/archive-software/github-upload/),
[DOI versioning](https://zenodo.org/help/versioning).

## PyPI and Bioconda

First publish a chosen version's exact source archive and wheel to PyPI when you
are ready (not performed here). Retain the exact source archive uploaded. The
future installation command will be `python -m pip install echidnaxomics`, followed
by `echidna --help`. Until publication, install from this checkout with `pip install -e .`.

Generate a local recipe draft from that exact sdist:

```bash
python packaging/render_bioconda_recipe.py dist/echidnaxomics-VERSION.tar.gz \
  --output-dir /tmp/echidnaxomics-bioconda-VERSION
```

Replace VERSION with the selected version. The script reads the version from
PKG-INFO, computes a real SHA256 and renders the noarch Python template. It does
not publish. Its URL becomes usable only after that exact archive is on PyPI;
a rebuilt archive can have a different hash even at the same version. Verify the
public URL/hash, copy the recipe to `recipes/echidnaxomics` in your fork of
`bioconda/bioconda-recipes`, run their build/lint checks, then open a PR when ready.
The recipe itself has not been built by conda/bioconda-utils locally.

After Bioconda accepts and builds it, installation will be:

```bash
conda create -n echidna -c conda-forge -c bioconda echidnaxomics
conda activate echidna
echidna --help
```

Install the required external tools separately as described in the installation
guide. They and their databases are not bundled in the base package.

References: [recipe guidelines](https://bioconda.github.io/contributor/guidelines.html),
[contribution workflow](https://bioconda.github.io/contributor/workflow.html).

## Wiki publication

The GitHub Wiki is a separate Git repository. Committing `wiki/Home.md` to the
code repository does not publish the Wiki. Follow [wiki activation and personal
commit/push commands](wiki.md); no publishing has been performed here.
