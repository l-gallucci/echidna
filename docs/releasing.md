# GitHub release 1.0.0: echidnaxomics

This is a preparation guide. No upload, DOI, release, or Bioconda submission has
been made by this preparation. The repository is still `l-gallucci/echidna`.
The selected software version is **1.0.0**, with intended tag **v1.0.0** and
GitHub Release title **Echidna 1.0.0**. [Release notes](releases/v1.0.0.md) are
ready to use. Creating a local tag or pushing it does not create a GitHub Release.

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
coauthors as additional `authors` entries. Version `1.0.0` and planned release
date `2026-11-04` are recorded. If publication occurs on a different date, update
`date-released` and the changelog date before tagging. The software version lives
in [`echidna/__init__.py`](../echidna/__init__.py); setup.py reads that value for
both the source archive and wheel.
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
no initial MMseqs classification, SSU or final-report and no chimera-based
exclusion. It can reuse a per-contig taxonomy table and resolve overlaps before
final QC, including recovery from dereplication clusters.

## Create the first GitHub Release

Commit and push the prepared changes yourself, then confirm that GitHub Actions
passes for that exact commit. Include all intended source, test and documentation
changes before tagging. A tag points to a committed revision and excludes any
uncommitted changes in the working tree.

```bash
git status --short
git add README.md CHANGELOG.md CITATION.cff MANIFEST.in mkdocs.yml \
  echidna/ tests/ docs/ wiki/ workflow/ examples/complete_workflow.sh \
  .github/workflows/ci.yml
git commit -m "Prepare Echidna 1.0.0 GitHub release"
git push origin main
```

After CI succeeds, confirm `git status --short` is empty and that HEAD is the
intended release commit, then tag and push only the new tag:

```bash
git tag -a v1.0.0 -m "Echidna 1.0.0"
git push origin v1.0.0
```

An existing local development tag `v2.0.0` was found during preparation. It was
not changed or deleted. The GitHub Releases list was empty on 2026-11-04;
the historical changelog is retained separately from the GitHub Releases list.
Do not push all local tags as part of this release.

Create a reviewable GitHub draft, using the exact tag already pushed:

```bash
gh release create v1.0.0 --repo l-gallucci/echidna --verify-tag --draft \
  --title "Echidna 1.0.0" --notes-file docs/releases/v1.0.0.md
```

Alternatively, open [New release](https://github.com/l-gallucci/echidna/releases/new),
select `v1.0.0`, use title `Echidna 1.0.0`, paste the notes and choose **Save draft**.
For the intended full release, leave **This is a pre-release** unchecked and
explicitly choose **Set as latest release** when publishing to label this as the
current release. Publication is a separate manual action; it has not been
performed here.

If you want the GitHub release to obtain a Zenodo DOI automatically, enable
and review the Zenodo integration **before publishing** the draft, as explained
below.

The wheel and source archive may be attached to the draft after building that
exact tagged revision. GitHub also provides source-code archives for the tag.
Neither action publishes to PyPI or Bioconda.

References: [GitHub release management](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository),
[GitHub CLI release creation](https://cli.github.com/manual/gh_release_create).

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
python packaging/render_bioconda_recipe.py dist/echidnaxomics-1.0.0.tar.gz \
  --output-dir /tmp/echidnaxomics-bioconda-1.0.0
```

The script reads the version from
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

## Documentation website

The documentation website is built with MkDocs and prepared for GitHub Pages.
It imports `wiki/Home.md` as a maintained source. Follow [website activation](wiki.md)
to select GitHub Actions as the Pages source once. Every push to `main` then
builds and publishes the site automatically; pull requests only validate it.
Documentation publication does not create a GitHub software release or issue a DOI.
