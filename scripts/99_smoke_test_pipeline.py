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
            '>sp|P1|PROT1\nMARRRRSTYKRRR\n'
            '>sp|P2|PROT2\nAKRSTYKKRRSTY\n'
            '>sp|P3|PROT3\nRRRSSSKKKYYYR\n'
        )
        base = pd.DataFrame([
            {'substrate_UniProtAC': 'P1', 'position': 3, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'unip'},
            {'substrate_UniProtAC': 'P1', 'position': 7, 'residue': 'S', 'ptm_group': 'Phosphorylation', 'is_arg_methyl': False, 'source': 'unip'},
            {'substrate_UniProtAC': 'P1', 'position': 10, 'residue': 'K', 'ptm_group': 'Ubiquitin/SUMO', 'is_arg_methyl': False, 'source': 'unip'},
            {'substrate_UniProtAC': 'P2', 'position': 4, 'residue': 'S', 'ptm_group': 'Phosphorylation', 'is_arg_methyl': False, 'source': 'unip'},
            {'substrate_UniProtAC': 'P2', 'position': 8, 'residue': 'K', 'ptm_group': 'Lys methylation', 'is_arg_methyl': False, 'source': 'unip'},
            {'substrate_UniProtAC': 'P3', 'position': 2, 'residue': 'R', 'ptm_group': 'Arg methylation', 'is_arg_methyl': True, 'source': 'iedb'},
            {'substrate_UniProtAC': 'P3', 'position': 5, 'residue': 'S', 'ptm_group': 'Phosphorylation', 'is_arg_methyl': False, 'source': 'unip'},
        ])
        integrated = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'corrected_position': 3, 'substrate_genename': 'GENE1', 'confidence_tier': 'Tier 1'},
            {'canonical_UniProtAC': 'P1', 'corrected_position': 5, 'substrate_genename': 'GENE1', 'confidence_tier': 'Tier 1'},
            {'canonical_UniProtAC': 'P3', 'corrected_position': 2, 'substrate_genename': 'GENE3', 'confidence_tier': 'Tier 2'},
            {'canonical_UniProtAC': 'P3', 'corrected_position': 4, 'substrate_genename': 'GENE3', 'confidence_tier': 'Tier 2'},
        ])
        condensate = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'is_condensate': True, 'role': 'driver', 'is_experimental': True, 'confidence_score': 0.95},
            {'canonical_UniProtAC': 'P3', 'is_condensate': True, 'role': 'member', 'is_experimental': False, 'confidence_score': 0.70},
        ])
        disorder = pd.DataFrame([
            {'canonical_UniProtAC': 'P1', 'fragment_start': 2, 'fragment_end': 6},
            {'canonical_UniProtAC': 'P3', 'fragment_start': 1, 'fragment_end': 5},
        ])

        write_tsv(base, tmp / 'base.tsv')
        write_tsv(integrated, tmp / 'integrated.tsv')
        write_tsv(condensate, tmp / 'condensate.tsv')
        write_tsv(disorder, tmp / 'disorder.tsv')

        assert nearest_distances_to_other([3, 5], [3, 5], exclude_self=True).tolist() == [2.0, 2.0]

        import subprocess
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
            'python', str(scripts / '12_cross_ptm_neighbor_enrichment.py'),
            '--base-master', str(tmp / 'base.tsv'),
            '--integrated-arg-sites', str(tmp / 'integrated.tsv'),
            '--canonical-fasta', str(fasta),
            '--disorder-intervals', str(tmp / 'disorder.tsv'),
            '--outdir', str(tmp / 'neighbors'),
            '--permutations', '10',
            '--max-distance', '10',
            '--progress-every', '0',
            '--checkpoints', '3', '5', '10',
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
            tmp / 'clustering' / 'nearest_neighbor_curve_vs_null.tsv',
            tmp / 'neighbors' / 'cross_ptm_neighbor_curves.tsv',
            tmp / 'condensate_out' / 'condensate_ptm_adjusted_logistic_models.tsv',
        ]
        for path in required:
            if not path.exists():
                raise FileNotFoundError(path)

        print({'status': 'ok', 'tmpdir': str(tmp)})


if __name__ == '__main__':
    main()
