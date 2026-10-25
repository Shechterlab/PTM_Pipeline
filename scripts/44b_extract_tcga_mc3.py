#!/usr/bin/env python3
"""Extract the TCGA somatic protein changes at arginines from the MC3 public MAF of the GDC.

MC3 (Ellrott et al. 2018, Cell Syst 6:271, doi:10.1016/j.cels.2018.03.002) called somatic mutations in the TCGA
exomes with seven callers; these are the calls of the TCGA PanCancer Atlas (Hoadley et al. 2018, Cell 173:291,
doi:10.1016/j.cell.2018.03.022). The GDC publishes the open-access calls as mc3.v0.2.8.PUBLIC.maf.gz
(https://gdc.cancer.gov/about-data/publications/mc3-2017), open-access data under the GDC's data access policies (NIH
Genomic Data Sharing policy), with no licence of their own.

A call is kept as the MC3 analyses kept them: FILTER PASS (no flag, and called by two callers or more), and in OV and
LAML also a call whose only flag is wga ("We restricted our analysis to PASS calls ... for OV and LAML we did not filter
mutations containing only 'wga' filter tags", Ellrott et al. 2018). Of those, the missense and nonsense changes at an
arginine of the MAF's transcript (HGVSp_Short p.R<position><residue or *>) are written; scripts/44_stage_depmap.py
keeps one only where the UniProt canonical sequence of the protein has an arginine at that position.

Which protein: UniProt (release 2026_03, --uniprot-ensembl and --uniprot-isoforms) names the Ensembl transcripts that
encode each reviewed entry's displayed (canonical) isoform. A call's change is read on such a transcript from the MAF's
all_effects (every transcript the variant falls in; mapped_by canonical_transcript). Only when none of the variant's
transcripts encodes a canonical isoform is the change of the MAF's own transcript (HGVSp_Short) written, without an
entry (mapped_by maf_transcript); scripts/44 then places it through the SWISSPROT entry the MAF names, or the gene symbol.

The project of each patient comes from the GDC cases endpoint (--cases: submitter_id, project.project_id,
project.name; fetched once and staged beside the output). A tumour is a sample (barcode to its sample type, e.g.
TCGA-02-0003-01): every tumour sample in the MAF counts as sequenced, with or without a kept call.

Outputs (in --outdir)
  tcga_mc3_arginine_changes_raw.tsv   one row per tumour sample, change and entry: project, sample, patient, gene, the
                                      entry (canonical transcript) or the UniProt entry name of the MAF transcript's
                                      protein, how it was mapped, change, transcript, FILTER, genomic change and read counts
  tcga_mc3_projects.tsv               each TCGA project: name, tumour samples in the MAF, those with a kept call, and
                                      the arginine changes written
  tcga_mc3_extract.json               inputs (MD5), the rule, and counts
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from common import md5_file

MAF_URL = 'https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc'
UNIPROT_ENSEMBL_URL = ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
                       '&fields=accession,id,xref_ensembl&format=tsv&compressed=true')
UNIPROT_ISOFORMS_URL = ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
                        '&fields=accession,cc_alternative_products&format=tsv&compressed=true')
CASES_URL = 'https://api.gdc.cancer.gov/cases (filter project.program.name = TCGA; fields submitter_id, project.project_id, project.name)'
CHANGE = re.compile(r'^p\.R(\d+)([ACDEFGHIKLMNPQSTVWY*])$')
CHANGE3 = re.compile(r'^p\.Arg(\d+)([A-Z][a-z]{2}|Ter)$')
THREE = {'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L', 'Met': 'M',
         'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'Ter': '*'}
KEEP_CONSEQUENCES = ('missense_variant', 'stop_gained')
KEEP_TYPES = {'Missense_Mutation', 'Nonsense_Mutation'}
WGA_PROJECTS = {'TCGA-OV', 'TCGA-LAML'}      # wga-only calls kept, as in Ellrott et al. 2018
FIELDS = ['project_id', 'sample_id', 'patient_id', 'hugo_symbol', 'canonical_UniProtAC', 'swissprot', 'mapped_by', 'protein_change',
          'mutation_type', 'transcript_id', 'filter', 'chromosome', 'start_position', 'reference_allele', 'tumor_allele',
          'tumor_alt_count', 'tumor_ref_count']


def canonical_transcripts(ensembl_tsv: str, isoforms_tsv: str) -> tuple[dict[str, set[str]], dict]:
    """Ensembl transcript (no version) -> the reviewed entries whose displayed isoform it encodes. A transcript UniProt
    tags with an isoform counts when that isoform is displayed; an untagged one only for an entry without isoforms."""
    displayed, counts = {}, Counter()
    with gzip.open(isoforms_tsv, 'rt', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            txt = r['Alternative products (isoforms)']
            if 'IsoId=' in txt:
                displayed[r['Entry']] = {x.strip() for d in re.findall(r'IsoId=([^;]+); Sequence=Displayed', txt) for x in d.split(',')}
    of = defaultdict(set)
    with gzip.open(ensembl_tsv, 'rt', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r['Entry']
            for item in filter(None, (x.strip() for x in r['Ensembl'].split(';'))):
                m = re.match(r'^(ENST\d+)(?:\.\d+)?(?: \[([A-Z0-9-]+)\])?$', item)
                if not m:
                    counts['xref_unparsed'] += 1
                    continue
                enst, iso = m.groups()
                if (iso and iso in displayed.get(acc, set())) or (not iso and acc not in displayed):
                    of[enst].add(acc)
                    counts['canonical_transcripts'] += 1
                else:
                    counts['other_isoform_transcripts'] += 1
    counts['entries_with_a_canonical_transcript'] = len({a for v in of.values() for a in v})
    return of, dict(counts)


def kept(filt: str, project: str) -> bool:
    return filt == 'PASS' or (project in WGA_PROJECTS and filt == 'wga')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--maf', required=True, help='mc3.v0.2.8.PUBLIC.maf.gz')
    ap.add_argument('--cases', required=True, help='GDC cases of TCGA (TSV with submitter_id, project.project_id, project.name)')
    ap.add_argument('--uniprot-ensembl', required=True, help='UniProt reviewed human: Entry, Entry Name, Ensembl (TSV, gzip)')
    ap.add_argument('--uniprot-isoforms', required=True, help='UniProt reviewed human: Entry, Alternative products (TSV, gzip)')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    project_of, project_name = {}, {}
    with open(args.cases, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            project_of[r['submitter_id']] = r['project.project_id']
            project_name[r['project.project_id']] = r['project.name']

    canon, xref_counts = canonical_transcripts(args.uniprot_ensembl, args.uniprot_isoforms)
    samples, samples_kept, written = defaultdict(set), defaultdict(set), Counter()
    counts = Counter()
    rows = []
    with gzip.open(args.maf, 'rt', encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        head = next(reader)
        ix = {c: k for k, c in enumerate(head)}
        for r in reader:
            counts['maf_rows'] += 1
            barcode = r[ix['Tumor_Sample_Barcode']]
            patient, sample = barcode[:12], barcode[:15]
            project = project_of.get(patient, '')
            if not project:
                counts['rows_patient_not_in_gdc_cases'] += 1
            samples[project].add(sample)
            filt = r[ix['FILTER']]
            if not kept(filt, project):
                counts['rows_filtered_out'] += 1
                continue
            samples_kept[project].add(sample)
            base = {'project_id': project, 'sample_id': sample, 'patient_id': patient, 'hugo_symbol': r[ix['Hugo_Symbol']],
                    'swissprot': r[ix['SWISSPROT']].replace('.', ''), 'filter': filt, 'chromosome': r[ix['Chromosome']],
                    'start_position': r[ix['Start_Position']], 'reference_allele': r[ix['Reference_Allele']],
                    'tumor_allele': r[ix['Tumor_Seq_Allele2']], 'tumor_alt_count': r[ix['t_alt_count']], 'tumor_ref_count': r[ix['t_ref_count']]}
            effects = [e.split(',') for e in r[ix['all_effects']].split(';') if e]
            on_canonical = [e for e in effects if len(e) >= 4 and e[3] in canon]
            if on_canonical:
                for e in on_canonical:
                    m = CHANGE3.match(e[2])
                    if not m or not any(c in e[1] for c in KEEP_CONSEQUENCES):
                        continue
                    for acc in sorted(canon[e[3]]):
                        rows.append(base | {'canonical_UniProtAC': acc, 'mapped_by': 'canonical_transcript', 'protein_change': f'R{m.group(1)}{THREE[m.group(2)]}',
                                            'mutation_type': e[1], 'transcript_id': e[3]})
                        counts['changes_at_R_on_a_canonical_transcript'] += 1
                        written[project] += 1
                continue
            if r[ix['Variant_Classification']] not in KEEP_TYPES:
                continue
            m = CHANGE.match(r[ix['HGVSp_Short']])
            if not m:
                continue
            counts['changes_at_R_on_the_maf_transcript_only'] += 1
            written[project] += 1
            rows.append(base | {'canonical_UniProtAC': '', 'mapped_by': 'maf_transcript', 'protein_change': r[ix['HGVSp_Short']][2:],
                                'mutation_type': r[ix['Variant_Classification']], 'transcript_id': r[ix['Transcript_ID']]})
    uniq = {tuple(x[f] for f in FIELDS): x for x in rows}   # one transcript of an entry per change (several can encode it)
    rows = sorted({(x['sample_id'], x['canonical_UniProtAC'], x['hugo_symbol'], x['protein_change'], x['chromosome'], x['start_position'], x['tumor_allele']): x
                   for x in sorted(uniq.values(), key=lambda x: x['transcript_id'])}.values(),
                  key=lambda x: (x['hugo_symbol'], x['canonical_UniProtAC'], int(CHANGE.match('p.' + x['protein_change']).group(1)), x['protein_change'],
                                 x['sample_id'], x['chromosome'], int(x['start_position']), x['tumor_allele']))
    with open(out / 'tcga_mc3_arginine_changes_raw.tsv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with open(out / 'tcga_mc3_projects.tsv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['project_id', 'name', 'tumour_samples', 'tumour_samples_with_kept_call', 'arginine_changes'])
        for p in sorted(x for x in samples if x):
            w.writerow([p, project_name.get(p, ''), len(samples[p]), len(samples_kept[p]), written[p]])
    meta = {'generated': dt.date.today().isoformat(), 'maf': {'file': Path(args.maf).name, 'path': str(Path(args.maf).resolve()), 'url': MAF_URL, 'md5': md5_file(args.maf)},
            'cases': {'file': Path(args.cases).name, 'url': CASES_URL, 'md5': md5_file(args.cases)},
            'uniprot': {'release': '2026_03', 'ensembl': {'file': Path(args.uniprot_ensembl).name, 'url': UNIPROT_ENSEMBL_URL, 'md5': md5_file(args.uniprot_ensembl)},
                        'isoforms': {'file': Path(args.uniprot_isoforms).name, 'url': UNIPROT_ISOFORMS_URL, 'md5': md5_file(args.uniprot_isoforms)},
                        'counts': xref_counts},
            'rule': 'FILTER PASS; in TCGA-OV and TCGA-LAML also FILTER wga alone (Ellrott et al. 2018); a missense or stop change at '
                    'an R of a transcript UniProt 2026_03 assigns to a reviewed entry\'s displayed isoform (all_effects); where the '
                    'variant falls in no such transcript, a Missense_Mutation or Nonsense_Mutation at an R of the MAF transcript '
                    '(HGVSp_Short)',
            'counts': dict(counts), 'tumour_samples': sum(len(v) for k, v in samples.items() if k),
            'projects': len([k for k in samples if k]), 'rows_written': len(rows)}
    (out / 'tcga_mc3_extract.json').write_text(json.dumps(meta, indent=1) + '\n')
    print(json.dumps(meta['counts']), meta['tumour_samples'], 'tumour samples,', len(rows), 'rows')


if __name__ == '__main__':
    main()
