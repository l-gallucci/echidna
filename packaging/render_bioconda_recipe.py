#!/usr/bin/env python3
"""Render a local Bioconda draft from the exact sdist intended for PyPI.

No network requests, uploads or submissions. The archive must later be uploaded
unchanged, otherwise its source hash in this recipe is no longer valid.
"""
import argparse
from email.parser import BytesParser
import hashlib
from pathlib import Path
import re
import tarfile


def render(archive, output_dir):
    with tarfile.open(archive, 'r:gz') as tar:
        metadata = [m for m in tar.getmembers()
                    if m.name.count('/') == 1 and m.name.endswith('/PKG-INFO') and m.isfile()]
        if len(metadata) != 1:
            raise ValueError('Expected one root PKG-INFO in the source archive')
        info = BytesParser().parsebytes(tar.extractfile(metadata[0]).read())
    if info['Name'] != 'echidnaxomics':
        raise ValueError('Expected an echidnaxomics source distribution')
    version = info['Version'] or ''
    if not re.fullmatch(r'[0-9]+(?:\.[0-9]+)*(?:(?:a|b|rc)[0-9]+)?(?:\.post[0-9]+)?(?:\.dev[0-9]+)?', version):
        raise ValueError('Unsupported version: ' + version)
    if archive.name != 'echidnaxomics-' + version + '.tar.gz':
        raise ValueError('Archive filename must match package name and version')
    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    template = (Path(__file__).parent / 'bioconda/meta.yaml.in').read_text()
    recipe = template.replace('@VERSION@', version).replace('@SHA256@', sha256)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / 'meta.yaml').write_text(recipe)
    print('Draft recipe:', output_dir / 'meta.yaml')
    print('Source SHA256:', sha256)
    print('Not submitted. Verify the public source URL after uploading this exact sdist.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    render(args.archive, args.output_dir)


if __name__ == '__main__':
    main()
