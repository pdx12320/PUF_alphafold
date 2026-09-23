"""Restore the supplied original archive after checking SHA256 and safe paths.

No downloads, training, Git operations or existing-file overwrites are performed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import stat
import zipfile
from pathlib import Path, PurePosixPath

EXPECTED_SHA256 = '09d5bcf49429ab849f716420d131b600d618f90874158c5257a862fbda8095c3'
PREFIX = 'c388_local_nonlocal_optimization'

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--out', type=Path, required=True,
                        help='New destination directory; must not exist')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Destination already exists; no files were overwritten')
    raw = args.archive.read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SHA256:
        raise ValueError('Archive SHA256 does not match the audited source')
    with zipfile.ZipFile(args.archive) as archive:
        infos = archive.infolist()
        if len(infos) != 94 or len({i.filename for i in infos}) != 94:
            raise ValueError('Unexpected archive member count or duplicate paths')
        for info in infos:
            path = PurePosixPath(info.filename)
            if (path.is_absolute() or '..' in path.parts or '\\' in info.filename
                    or path.parts[0] != PREFIX or stat.S_ISLNK(info.external_attr >> 16)):
                raise ValueError(f'Unsafe member: {info.filename}')
        if archive.testzip() is not None:
            raise ValueError('ZIP CRC check failed')
        manifest = json.loads(archive.read(f'{PREFIX}/artifact_sha256.json'))
        if len(manifest) != 93:
            raise ValueError('Unexpected source manifest size')
        for relative, digest in manifest.items():
            if hashlib.sha256(archive.read(f'{PREFIX}/{relative}')).hexdigest() != digest:
                raise ValueError(f'Checksum mismatch: {relative}')
        args.out.mkdir(parents=True, exist_ok=False)
        archive.extractall(args.out)
    print(f'Restored 94 original files to {args.out / PREFIX}')

if __name__ == '__main__':
    main()
