from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd

from clustering_utils import nearest_distances_to_other


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep='\t', index=False)


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        scripts = Path(__file__).resolve().parent
        repo = scripts.parent

        fasta = tmp / 'mini.fasta'
        fasta.write_text(
            '>sp|P1|FUS_HUMAN\nMGRGGGRGGPRGGGFGGDRGGYGGSRGGGPA\n'
            '>sp|P2|DDX3X_HUMAN\nMARRGTPQDSAARRGGKDVVTRRGAAPG\n'
            '>sp|P3|RBMX_HUMAN\nMRGGRSPPPGGRGGRDDDRGGRAAAPGG\n'
            '>sp|P4|KMT2D_HUMAN\nMSTAAARRKSTEEEDDRKRRPTAAAGG\n'
        )
        base = pd.DataFrame([
            {'substrate_UniProtAC': 'P1', 'position': 2, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'unip', 'protein_name': 'FUS RNA-binding protein', 'UniProtID': 'FUS_HUMAN', 'substrate_genename': 'FUS'},
            {'substrate_UniProtAC': 'P1', 'position': 6, 'residue': 'S', 'ptm_group': 'Phosphorylation', 'is_arg_methyl': False, 'source': 'unip', 'protein_name': 'FUS RNA-binding protein', 'UniProtID': 'FUS_HUMAN', 'substrate_genename': 'FUS'},
            {'substrate_UniProtAC': 'P2', 'position': 4, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'iedb', 'protein_name': 'ATP-dependent RNA helicase DDX3X', 'UniProtID': 'DDX3X_HUMAN', 'substrate_genename': 'DDX3X'},
            {'substrate_UniProtAC': 'P2', 'position': 9, 'residue': 'S', 'ptm_group': 'Phosphorylation', 'is_arg_methyl': False, 'source': 'unip', 'protein_name': 'ATP-dependent RNA helicase DDX3X', 'UniProtID': 'DDX3X_HUMAN', 'substrate_genename': 'DDX3X'},
            {'substrate_UniProtAC': 'P3', 'position': 3, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'unip', 'protein_name': 'RNA-binding motif protein X-linked spliceosome factor', 'UniProtID': 'RBMX_HUMAN', 'substrate_genename': 'RBMX'},
            {'substrate_UniProtAC': 'P3', 'position': 12, 'residue': 'K', 'ptm_group': 'Lys methylation', 'is_arg_methyl': False, 'source': 'unip', 'protein_name': 'RNA-binding motif protein X-linked spliceosome factor', 'UniProtID': 'RBMX_HUMAN', 'substrate_genename': 'RBMX'},
            {'substrate_UniProtAC': 'P4', 'position': 8, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'unip', 'protein_name': 'Histone-lysine N-methyltransferase KMT2D chromatin regulator', 'UniProtID': 'KMT2D_HUMAN', 'substrate_genename': 'KMT2D'},
            {'substrate_UniProtAC': 'P4', 'position': 15, 'residue': 'K', 'ptm_group': 'Ubiquitin/SUMO', 'is_arg_methyl': False, 'source': 'unip', 'protein_name': 'Histone-lysine N-methyltransferase KMT2D chromatin regulator', 'UniProtID': 'KMT2D_HUMAN', 'substrate_genename': 'KMT2D'},
        ])
        integrated = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'corrected_position': 2, 'residue': 'R', 'substrate_genename': 'FUS', 'site_key': 'P1:R2'},
            {'canonical_UniProtAC': 'P1', 'corrected_position': 5, 'residue': 'R', 'substrate_genename': 'FUS', 'site_key': 'P1:R5'},
            {'canonical_UniProtAC': 'P2', 'corrected_position': 4, 'residue': 'R', 'substrate_genename': 'DDX3X', 'site_key': 'P2:R4'},
            {'canonical_UniProtAC': 'P3', 'corrected_position': 3, 'residue': 'R', 'substrate_genename': 'RBMX', 'site_key': 'P3:R3'},
            {'canonical_UniProtAC': 'P4', 'corrected_position': 8, 'residue': 'R', 'substrate_genename': 'KMT2D', 'site_key': 'P4:R8'},
        ])
        condensate = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'is_condensate': True, 'role': 'driver', 'is_experimental': True, 'confidence_score': 0.95},
            {'canonical_UniProtAC': 'P2', 'is_condensate': True, 'role': 'member', 'is_experimental': True, 'confidence_score': 0.80},
            {'canonical_UniProtAC': 'P3', 'is_condensate': True, 'role': 'member', 'is_experimental': False, 'confidence_score': 0.70},
        ])
        disorder = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'fragment_start': 1, 'fragment_end': 20},
            {'canonical_UniProtAC': 'P2', 'fragment_start': 1, 'fragment_end': 12},
            {'canonical_UniProtAC': 'P3', 'fragment_start': 1, 'fragment_end': 18},
        ])
        interpro = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'interpro_accession': 'IPR0001', 'interpro_name': 'RGG low complexity RNA-binding region', 'interpro_type': 'domain', 'member_database_accessions': 'RGG', 'fragment_start': 1, 'fragment_end': 10},
            {'canonical_UniProtAC': 'P2', 'interpro_accession': 'IPR0002', 'interpro_name': 'DEAD-box helicase', 'interpro_type': 'domain', 'member_database_accessions': 'DDX', 'fragment_start': 1, 'fragment_end': 12},
            {'canonical_UniProtAC': 'P3', 'interpro_accession': 'IPR0003', 'interpro_name': 'RNA recognition motif', 'interpro_type': 'domain', 'member_database_accessions': 'RRM', 'fragment_start': 1, 'fragment_end': 8},
            {'canonical_UniProtAC': 'P4', 'interpro_accession': 'IPR0004', 'interpro_name': 'SET chromatin methyltransferase domain', 'interpro_type': 'domain', 'member_database_accessions': 'SET', 'fragment_start': 10, 'fragment_end': 22},
        ])

        write_tsv(base, tmp / 'base.tsv')
        write_tsv(integrated, tmp / 'integrated.tsv')
        write_tsv(condensate, tmp / 'condensate.tsv')
        write_tsv(disorder, tmp / 'disorder.tsv')
        write_tsv(interpro, tmp / 'interpro.tsv')
        (tmp / 'empty_interpro.tsv').write_text('')

        assert nearest_distances_to_other([3, 5], [3, 5], exclude_self=True).tolist() == [2.0, 2.0]

        import subprocess
        subprocess.run([
            'python', str(scripts / '03_functional_class_enrichment.py'),
            '--base-master', str(tmp / 'base.tsv'),
            '--integrated-sites', str(tmp / 'integrated.tsv'),
            '--outdir', str(tmp / 'functional'),
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '05_annotate_domain_context.py'),
            '--sites', str(tmp / 'integrated.tsv'),
            '--interpro-intervals', str(tmp / 'interpro.tsv'),
            '--outdir', str(tmp / 'domain_annot'),
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '05_annotate_domain_context.py'),
            '--sites', str(tmp / 'integrated.tsv'),
            '--interpro-intervals', str(tmp / 'empty_interpro.tsv'),
            '--outdir', str(tmp / 'domain_empty'),
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '07_summarize_domain_enrichment.py'),
            '--annotated-sites', str(tmp / 'domain_annot' / 'sites_with_domain_context.tsv'),
            '--interpro-intervals', str(tmp / 'interpro.tsv'),
            '--canonical-fasta', str(fasta),
            '--outdir', str(tmp / 'domain_summary'),
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '08_motif_and_arg_odds.py'),
            '--integrated-sites', str(tmp / 'integrated.tsv'),
            '--canonical-fasta', str(fasta),
            '--disorder-intervals', str(tmp / 'disorder.tsv'),
            '--outdir', str(tmp / 'motif'),
            '--min-kmer-count', '1',
            '--min-arg-count', '1',
            '--min-site-count', '1',
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '10_methyl_arg_clustering.py'),
            '--integrated-sites', str(tmp / 'integrated.tsv'),
            '--canonical-fasta', str(fasta),
            '--outdir', str(tmp / 'clustering'),
            '--permutations', '10',
            '--max-distance', '10',
            '--progress-every', '0',
        ], check=True, cwd=repo)
        subprocess.run([
            'python', str(scripts / '13_condensate_ptm_enrichment.py'),
            '--base-master', str(tmp / 'base.tsv'),
            '--integrated-arg-sites', str(tmp / 'integrated.tsv'),
            '--canonical-fasta', str(fasta),
            '--condensate-proteins', str(tmp / 'condensate.tsv'),
            '--disorder-intervals', str(tmp / 'disorder.tsv'),
            '--outdir', str(tmp / 'condensate_out'),
        ], check=True, cwd=repo)

        required = [
            tmp / 'functional' / 'arg_methyl_functional_multilabel_enrichment.tsv',
            tmp / 'domain_annot' / 'sites_with_domain_context.tsv',
            tmp / 'domain_empty' / 'sites_with_domain_context.tsv',
            tmp / 'domain_summary' / 'arg_methyl_domain_edge_proximity_curve.pdf',
            tmp / 'motif' / 'exclusive_motif_family_enrichment_vs_matched_nonmethyl_arginines.tsv',
            tmp / 'motif' / 'motif_feature_model.tsv',
            tmp / 'clustering' / 'nearest_neighbor_curve_vs_null.tsv',
            tmp / 'condensate_out' / 'condensate_ptm_adjusted_logistic_models.tsv',
        ]
        for path in required:
            if not path.exists():
                raise FileNotFoundError(path)

        print({'status': 'ok', 'tmpdir': str(tmp)})


if __name__ == '__main__':
    main()
