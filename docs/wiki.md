# Documentation website

The documentation website is prepared for **GitHub Pages** at
[https://l-gallucci.github.io/echidna/](https://l-gallucci.github.io/echidna/).
It becomes available after you enable Pages and run the publication workflow.
No website has been published as part of this preparation.

## Activate GitHub Pages

After reviewing, committing and pushing these files yourself:

1. Open the repository's **Settings → Pages**.
2. Under **Build and deployment → Source**, select **GitHub Actions**.
3. Open **Actions → Documentation website → Run workflow**.
4. Choose the `main` branch, enable **Publish the built documentation to GitHub Pages**,
   and click **Run workflow**.
5. When the `deploy` job succeeds, visit the URL reported by that job.

Pushes and pull requests build and validate the documentation. Publication is
manual: when you want to update the live site, run the same workflow with the
publish option enabled. A run without that option only validates the site.
Publishing documentation does not create a software release or DOI.

## Edit the content

| Content | Maintained source |
|---|---|
| Website home | `docs/index.md` |
| Installation, workflows, commands, references | Corresponding Markdown files in `docs/` |
| Interpretation guide (Italian) | `wiki/Home.md` |
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
