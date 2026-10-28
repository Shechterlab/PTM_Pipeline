import sys, json
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import lib
repo, src, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ctx = lib.setup(repo)
t = pd.read_csv(src / 'data/labels/cluster_annotation_Dec6.csv', dtype=str, keep_default_na=False)
t['acc'] = [lib.resolve(ctx, ids.split(';'), g) for ids, g in zip(t['Majority protein IDs'], t['Gene_name_canonical'])]
lost = (t['acc'] == '').sum()
t = t[t['acc'] != ''].drop_duplicates('acc')
loc = t[t['cluster_annotation'] != 'unclassified']
print(f'groups {len(t) + lost}; on a reviewed entry {len(t)}; unmapped {lost}; classified {len(loc)}')
m = ctx['meth']
print(f'methylated: detected {len(set(t.acc) & m)} ({len(set(t.acc) & m)/len(t):.1%}); classified {len(set(loc.acc) & m)} ({len(set(loc.acc) & m)/len(loc):.1%})')
sets = {c: set(g['acc']) for c, g in loc.groupby('cluster_annotation')}
res = lib.tests(ctx, set(loc.acc), sets)
lib.show(res, 'organelle IP (Hein 2025), classified proteins')
res.to_csv(out / 'oip_pilot.tsv', sep='\t', index=False)
json.dump(dict(zip(t.acc, t.cluster_annotation)), open(out / 'oip_protein.json', 'w'))
