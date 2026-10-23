"""Over-representation tests with published tools: clusterProfiler::enricher and DOSE::enrichDGN, run by
scripts/enrich_sets.R (Rscript on PATH, from the pipeline's environment).

    from enrichment import Job, run
    results = run([Job('panther_group', targets, universe, sets=term2gene, whole_universe=True)], workdir)

Each result is clusterProfiler's table for every tested set (one-sided hypergeometric pvalue, BH p.adjust within the
job, FoldEnrichment), with the counts of its ratios as integers: k target identifiers in the set, n targets in the
background, M background identifiers in the set, N background, and the 2x2 table a = k, b = n - k, c = M - k,
d = N - M - n + k.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

RUNNER = Path(__file__).with_name('enrich_sets.R')


@dataclass
class Job:
    job: str
    target: object                      # iterable of identifiers
    universe: object                    # iterable of identifiers
    sets: pd.DataFrame | None = None    # term, gene (TERM2GENE)
    names: pd.DataFrame | None = None   # term, name (TERM2NAME)
    whole_universe: bool = False
    min_size: int = 1
    max_size: int = 10 ** 6
    tool: str = 'enricher'
    extra: dict = field(default_factory=dict)


def inputs_key(files: dict[str, str], row: dict, runner: bytes) -> str:
    """What a job's results depend on: its input files' text, its settings and the R script."""
    h = hashlib.sha256(runner)
    for name in sorted(files):
        h.update(f'\0{name}\0{files[name]}'.encode())
    h.update(json.dumps({k: str(v) for k, v in row.items() if k not in ('target', 'universe', 'sets', 'names')}, sort_keys=True).encode())
    return h.hexdigest()


def run(jobs: list[Job], workdir: Path) -> dict[str, pd.DataFrame]:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / 'results'
    runner = RUNNER.read_bytes()
    rows, keys = [], {}
    for j in jobs:
        d = workdir / j.job
        d.mkdir(exist_ok=True)
        files = {'target.txt': '\n'.join(sorted({str(x) for x in j.target})) + '\n',
                 'universe.txt': '\n'.join(sorted({str(x) for x in j.universe})) + '\n'}
        if j.sets is not None:
            files['sets.tsv'] = j.sets.iloc[:, :2].set_axis(['term', 'gene'], axis=1).drop_duplicates().to_csv(sep='\t', index=False)
        if j.names is not None:
            files['names.tsv'] = j.names.iloc[:, :2].set_axis(['term', 'name'], axis=1).drop_duplicates().to_csv(sep='\t', index=False)
        for name, text in files.items():
            (d / name).write_text(text)
        row = {'job': j.job, 'tool': j.tool, 'target': str(d / 'target.txt'), 'universe': str(d / 'universe.txt'),
               'sets': str(d / 'sets.tsv') if 'sets.tsv' in files else '', 'names': str(d / 'names.tsv') if 'names.tsv' in files else '',
               'whole_universe': str(j.whole_universe).lower(), 'min_size': j.min_size, 'max_size': j.max_size}
        rows.append(row)
        keys[j.job] = inputs_key(files, row, runner)
    # one jobs table and log per call, named after its jobs: several calls can share a workdir
    tag = '_'.join(j.job for j in jobs)[:120]
    jobs_file, log = workdir / f'jobs_{tag}.tsv', workdir / f'enrich_sets_{tag}.log'
    jobs_file.write_text(pd.DataFrame(rows).to_csv(sep='\t', index=False))
    # A job's results are reused when the stamp R left beside them (written only after it succeeded) names the same
    # inputs, settings and script; otherwise R runs again.
    stamp = lambda job: out / f'{job}.inputs.sha256'  # noqa: E731
    if not ((out / 'versions.tsv').exists() and all((out / f'{job}.tsv').exists() and stamp(job).exists()
                                                    and stamp(job).read_text().strip() == key for job, key in keys.items())):
        rscript = shutil.which('Rscript')
        if not rscript:
            raise SystemExit('Rscript is not on PATH: the enrichment tests need R with clusterProfiler and DOSE (environment file)')
        for job in keys:
            stamp(job).unlink(missing_ok=True)
        proc = subprocess.run([rscript, str(RUNNER), str(jobs_file), str(out)], capture_output=True, text=True)
        log.write_text(proc.stdout + proc.stderr)
        if proc.returncode:
            raise SystemExit(f'scripts/enrich_sets.R failed (log: {log}):\n{proc.stderr[-2000:]}')
        for job, key in keys.items():
            stamp(job).write_text(key + '\n')
    results = {}
    for j in jobs:
        t = pd.read_csv(out / f'{j.job}.tsv', sep='\t', dtype={'ID': str, 'Description': str, 'geneID': str}, keep_default_na=False, na_values=[''])
        if len(t):
            t[['k', 'n']] = t['GeneRatio'].str.split('/', expand=True).astype(int)
            t[['M', 'N']] = t['BgRatio'].str.split('/', expand=True).astype(int)
            t['a'], t['b'], t['c'] = t['k'], t['n'] - t['k'], t['M'] - t['k']
            t['d'] = t['N'] - t['M'] - t['b']
        results[j.job] = t
    return results


def versions(workdir: Path) -> dict[str, str]:
    t = pd.read_csv(Path(workdir) / 'results' / 'versions.tsv', sep='\t', dtype=str)
    return dict(zip(t['package'], t['version']))
