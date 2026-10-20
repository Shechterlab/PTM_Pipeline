#!/usr/bin/env python3
"""Create a Zenodo draft deposit, upload files to it and set its metadata. It never publishes: check the draft on
zenodo.org and press Publish there.

    export ZENODO_TOKEN=...      # zenodo.org > Account > Applications > Personal access token (scopes deposit:write, deposit:actions)
    python zenodo_upload.py --metadata zenodo_public.json \
        README.md MD5SUMS SHA256SUMS \
        ptm_pipeline_inputs_tables_2026-10-20.tar.gz \
        ptm_pipeline_interactions_raw_cache_2026-10-16.tar.gz \
        ptm_pipeline_layers_raw_downloads_2026-10-18.tar.gz

--sandbox uses sandbox.zenodo.org (a separate account and token) for a trial run. --deposit ID adds files to an
existing draft instead of creating one. Each upload is checked against the local MD5 that Zenodo returns.
Standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path


def call(method: str, url: str, token: str, body=None, headers=None, length: int | None = None):
    h = {'Authorization': f'Bearer {token}', **(headers or {})}
    if length is not None:
        h['Content-Length'] = str(length)
    req = urllib.request.Request(url, data=body, method=method, headers=h)
    try:
        with urllib.request.urlopen(req, timeout=3600) as r:
            text = r.read().decode()
            return json.loads(text) if text else {}
    except urllib.error.HTTPError as e:
        sys.exit(f'{method} {url}: HTTP {e.code}\n{e.read().decode(errors="replace")[:2000]}')


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--metadata', required=True, help='JSON file with {"metadata": {...}} (zenodo_public.json or zenodo_restricted.json)')
    ap.add_argument('--sandbox', action='store_true', help='use sandbox.zenodo.org')
    ap.add_argument('--deposit', help='add to this existing draft deposit id instead of creating a new one')
    ap.add_argument('--token-env', default='ZENODO_TOKEN')
    ap.add_argument('files', nargs='+', type=Path)
    args = ap.parse_args()

    token = os.environ.get(args.token_env)
    if not token:
        sys.exit(f'set {args.token_env} to a Zenodo personal access token')
    missing = [str(f) for f in args.files if not f.is_file()]
    if missing:
        sys.exit(f'missing files: {missing}')
    api = 'https://sandbox.zenodo.org/api' if args.sandbox else 'https://zenodo.org/api'
    metadata = json.loads(Path(args.metadata).read_text())

    if args.deposit:
        dep = call('GET', f'{api}/deposit/depositions/{args.deposit}', token)
    else:
        dep = call('POST', f'{api}/deposit/depositions', token, b'{}', {'Content-Type': 'application/json'})
        print(f'created draft deposit {dep["id"]}')
    bucket = dep['links']['bucket']

    for f in args.files:
        size = f.stat().st_size
        print(f'uploading {f.name} ({size / 1e6:.1f} MB) ...', flush=True)
        with open(f, 'rb') as fh:
            res = call('PUT', f'{bucket}/{f.name}', token, fh, {'Content-Type': 'application/octet-stream'}, length=size)
        remote = (res.get('checksum') or '').replace('md5:', '')
        local = md5(f)
        if remote != local:
            sys.exit(f'{f.name}: Zenodo reports MD5 {remote}, local file is {local}')
        print(f'  ok, md5 {local}')

    call('PUT', f'{api}/deposit/depositions/{dep["id"]}', token, json.dumps(metadata).encode(), {'Content-Type': 'application/json'})
    site = 'https://sandbox.zenodo.org' if args.sandbox else 'https://zenodo.org'
    print(f'metadata set. Review and publish the draft at {site}/deposit/{dep["id"]}')


if __name__ == '__main__':
    main()
