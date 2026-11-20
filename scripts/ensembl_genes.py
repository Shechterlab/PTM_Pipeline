"""Gene names for the variant tables: Ensembl transcripts to genes, and every gene to its current HGNC approved symbol.

The four sources name genes in different years' terms (the MC3 MAF in 2017 symbols on GENCODE v19 transcripts, DepMap and
ClinVar in their own, gnomAD v4.1 by Ensembl 105), so each row's gene is resolved to the HGNC approved symbol through a stable
identifier where the source gives one (Ensembl gene, NCBI Gene), else through the symbol it gives (approved, then a unique
previous symbol, then a unique alias). HGNC data are CC0 (https://genenames.org/about/license).
"""
from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict

GTF_105 = 'https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz'
HGNC_URL = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'


def transcript_genes(gtf: str) -> dict[str, tuple[str, str]]:
    """Ensembl transcript (no version) -> (Ensembl gene, gene name) in an Ensembl GTF."""
    out: dict[str, tuple[str, str]] = {}
    with gzip.open(gtf, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if len(f) < 9 or f[2] != 'transcript':
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8])
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            name = re.search(r'gene_name "([^"]+)"', f[8])
            if tid:
                out[tid.group(1).split('.')[0]] = (gid.group(1).split('.')[0] if gid else '', name.group(1) if name else '')
    return out


class HGNC:
    """The HGNC complete set (approved genes): symbols by Ensembl gene, NCBI Gene and by any symbol a gene has had."""

    def __init__(self, path: str):
        self.by_ensg, self.by_entrez, self.approved = {}, {}, set()
        prev, alias = defaultdict(set), defaultdict(set)
        with open(path, encoding='utf-8', newline='') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['status'] != 'Approved':
                    continue
                s = r['symbol']
                self.approved.add(s)
                if r['ensembl_gene_id']:
                    self.by_ensg[r['ensembl_gene_id']] = s
                if r['entrez_id']:
                    self.by_entrez[r['entrez_id']] = s
                for p in filter(None, r['prev_symbol'].strip('"').split('|')):
                    prev[p].add(s)
                for a in filter(None, r['alias_symbol'].strip('"').split('|')):
                    alias[a].add(s)
        self.prev = {k: next(iter(v)) for k, v in prev.items() if len(v) == 1}
        self.alias = {k: next(iter(v)) for k, v in alias.items() if len(v) == 1}

    def from_symbol(self, name: str) -> str:
        out = []
        for part in filter(None, (x.strip() for x in re.split(r'[;,]', name or ''))):
            s = part if part in self.approved else self.prev.get(part) or self.alias.get(part) or part
            if s not in out:
                out.append(s)
        return ';'.join(out)

    def symbol(self, ensg: str = '', entrez: str = '', name: str = '') -> str:
        ensg = (ensg or '').split('.')[0]
        if ensg in self.by_ensg:
            return self.by_ensg[ensg]
        if entrez in self.by_entrez:
            return self.by_entrez[entrez]
        return self.from_symbol(name)
