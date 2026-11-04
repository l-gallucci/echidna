# Documentation website

The documentation website is prepared for **GitHub Pages** at
[https://l-gallucci.github.io/echidna/](https://l-gallucci.github.io/echidna/).
The workflow builds and publishes it automatically after every push to `main`.

## Activate GitHub Pages

After reviewing, committing and pushing these files yourself:

1. Open the repository's **Settings → Pages**.
2. Under **Build and deployment → Source**, select **GitHub Actions**.
3. Push a commit to `main`. Alternatively, open **Actions → Documentation website
   → Run workflow**, select `main`, and run it to redeploy the current version.
4. When the `deploy` job succeeds, visit the URL reported by that job.

Every push to `main` builds, validates and deploys the site, including pushes
that only change code. Pull requests build and validate without deploying.
Manual runs deploy only when `main` is selected. A failed build prevents
deployment; the previously deployed site remains available. No publish checkbox
is required. If Pages is already configured with GitHub Actions, no further
activation is needed.
Publishing documentation does not create a software release or DOI.

## Edit the content

| Content | Maintained source |
|---|---|
| Website home | `docs/index.md` |
| Installation, workflows, commands, references | Corresponding Markdown files in `docs/` |
| Interpretation guide (English) | `wiki/Home.md` |
| Snakemake guide | `workflow/README.md` |
| Changelog | `CHANGELOG.md` |
| Navigation and theme | `mkdocs.yml` |

The build imports the wiki guide, Snakemake guide and changelog directly from
their maintained sources. You do not need to edit a second copy. The audit
records remain in the repository and are excluded from the public website.

## Local preview

Use Python 3.12 in a separate documentation environment; these packages are not required by
the `echidna` command:

```bash
python3.12 -m venv .venv-docs
source .venv-docs/bin/activate
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

Open the local URL printed by MkDocs. To check the production build:

```bash
python -m mkdocs build --strict
```

Generated output goes to the ignored `site/` directory. Commit the Markdown,
configuration, hook and workflow files, rather than generated HTML.

## GitHub Wiki

The repository's integrated GitHub Wiki is a separate, optional service. It is
not required for this website. `wiki/Home.md` is used here as a Markdown source;
there is no need to clone or push a separate `.wiki.git` repository.

Official references: [custom GitHub Pages workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages),
[MkDocs](https://www.mkdocs.org/).
