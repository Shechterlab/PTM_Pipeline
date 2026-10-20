# Human methylarginine sites, v1.0.0 (published)

The site union published with DeAngelo, Silverstein & Shechter, *Protein arginine methylation rewires molecular
interactions*, Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011: 13,719 sites on 5,090
accessions, from 28,960 source records (28,134 placed).

These two tables let later releases say exactly what changed. Every value is the published one; only the column
names follow later releases.

| File | One row per |
|---|---|
| `human_methylarginine_sites_v1.0.0.tsv` | published site: accession, gene, position, methylation states, source tables, number of records, and the disorder and domain classes the published analyses used |
| `supporting/source_records_v1.0.0.tsv` | source record: what the study reported, its evidence sequence, and where the published run placed it (`placed_site`); `record_id` matches later releases |
| `MD5SUMS` | checksums (`md5sum -c MD5SUMS`) |

**Where the full published results are.** The analysis outputs, figures and website data of v1.0.0 are in the
repository at git tag `v1.0.0` (folder `PTM_results/`) and in the `v1.0.0_archive` branch, and the frozen website is
in `docs/v1.0.0/`. `scripts/build_release_v1_0_0.py` wrote these tables from that folder.

**Reference sequences.** Positions refer to the UniProtKB reviewed human sequences staged on 2026-04-07, whose
UniProt release was not recorded (it predates 2026_03). Some published sites sit on unreviewed or obsolete
accessions; v1.1.0 corrects these (see `../v1.1.0/CHANGES.md`).
