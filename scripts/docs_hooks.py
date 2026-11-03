"""Build the website from the maintained docs and wiki, without duplicate sources."""
import posixpath
from pathlib import Path
import re
from urllib.parse import quote, unquote, urlsplit, urlunsplit

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
GENERATED = {
    'guide.md': ROOT / 'wiki/Home.md',
    'snakemake.md': ROOT / 'workflow/README.md',
    'changelog.md': ROOT / 'CHANGELOG.md',
}
SOURCE_TO_PAGE = {source.resolve(): uri for uri, source in GENERATED.items()}
SOURCE_TO_PAGE[ROOT / 'README.md'] = 'index.md'
SOURCE_TO_PAGE[ROOT / 'INSTALL.md'] = 'installation.md'
LINK = re.compile(r'(!?\[[^\]]*\]\()([^\s)]+)([^)]*\))')


def on_files(files, config):
    for uri, source in GENERATED.items():
        files.append(File.generated(config, uri, content=source.read_text(encoding='utf-8')))
    files.append(File.generated(config, 'img/logo.jpg', content=(ROOT / 'assets/LOGO2.jpg').read_bytes()))
    return files


def on_page_markdown(markdown, page, config, files):
    """Resolve repository links correctly for the published project subdirectory."""
    source = GENERATED.get(page.file.src_uri, DOCS / page.file.src_uri)

    def rewrite(match):
        target = urlsplit(match.group(2))
        if target.scheme or target.netloc or not target.path or target.path.startswith('/'):
            return match.group(0)
        candidate = (source.parent / unquote(target.path)).resolve()
        destination = SOURCE_TO_PAGE.get(candidate)
        if destination is None and candidate.is_relative_to(DOCS):
            destination = candidate.relative_to(DOCS).as_posix()
        if destination is not None:
            relative = posixpath.relpath(destination, posixpath.dirname(page.file.src_uri) or '.')
            url = urlunsplit(('', '', relative, target.query, target.fragment))
        elif candidate.is_relative_to(ROOT) and candidate.exists():
            url = 'https://github.com/l-gallucci/echidna/blob/main/' + quote(candidate.relative_to(ROOT).as_posix())
            if target.fragment:
                url += '#' + target.fragment
        else:
            return match.group(0)
        return match.group(1) + url + match.group(3)

    # Leave copy/paste shell examples unchanged.
    parts = re.split(r'(^```[^\n]*\n.*?^```\s*$)', markdown, flags=re.MULTILINE | re.DOTALL)
    return ''.join(part if part.startswith('```') else LINK.sub(rewrite, part) for part in parts)
