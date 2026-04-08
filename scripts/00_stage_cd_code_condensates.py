from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import pandas as pd

from common import save_table


def load_condensate_module():
    script_path = Path(__file__).with_name('13_condensate_ptm_enrichment.py')
    spec = importlib.util.spec_from_file_location('condensate_ptm_enrichment', script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'could not load module from {script_path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--cdcode-proteins', required=True)
    ap.add_argument('--cdcode-protein-map', default='')
    ap.add_argument('--cdcode-condensates', default='')
    ap.add_argument('--out', default='data/context/cd_code_condensate_proteins.tsv')
    args = ap.parse_args()

    module = load_condensate_module()
    staged = module.build_condensate_table_from_cd_code_raw(
        proteins_path=args.cdcode_proteins,
        protein_map_path=args.cdcode_protein_map,
        condensates_path=args.cdcode_condensates,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    save_table(staged, out)
    summary = pd.DataFrame([{
        'proteins_total': int(len(staged)),
        'condensate_positive': int(staged['is_condensate'].fillna(False).sum()),
        'output': str(out),
    }])
    save_table(summary, out.with_name(out.stem + '_summary.tsv'))
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()
