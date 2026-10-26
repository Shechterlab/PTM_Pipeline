# Adding a data source

There are two kinds of source:

- **A study or database of methylarginine sites** adds sites. It is staged as one table of records by
  `scripts/00_stage_study_sites.py` into `inputs/sources/`, listed in `STUDIES` in `scripts/release_codes.py`, and read
  by `hpc/build_release.sh` (`STUDIES`). Every record then goes through the same placement and site tests as the rest.
- **An annotation layer** describes sites or proteins that are already there, and adds no site. Examples: the CARM1
  knockout response, Human Protein Atlas locations, DepMap and TCGA. The rest of this file is about layers.

A layer goes through four steps: stage it, register it, show it on the website, check it.

## 1. Stage it

Write `scripts/NN_stage_<name>.py`. It reads the downloads and writes tables into `inputs/layers/`. `scripts/44_stage_depmap.py`
is the model to copy. Keep to these rules:

- **Keys.** Proteins are keyed by `canonical_UniProtAC`, the reviewed accession of the reference proteome
  (`inputs/reference/uniprot_human_reviewed_2026_03.fasta`). Residue data also carry `position` on that sequence.
- **Residues.** Keep a residue record only where the canonical sequence has that residue: `common.residue_is`. Count
  what is dropped and why in the summary JSON.
- **Genes.** To map gene-level data (Entrez genes) to entries, use `common.entrez_of_entries`. It maps them the way the
  disease tests do.
- **Provenance.** Write a summary JSON with counts, versions and the MD5 of every input. Each download goes into
  `inputs/layers/raw_files.tsv` (URL, size, MD5, SHA-256, date) with `common.upsert_tsv_rows`. Then refresh that file's
  own MANIFEST row with `common.refresh_manifest_row`.
- **Manifest.** Each staged table gets a row in `inputs/MANIFEST.tsv` (`common.update_manifest`). The row gives its
  content, source, version, citation, licence, how it was made, and `redistribution`:
  - `git` if the licence allows the table in the repository;
  - `zenodo` for raw downloads kept in the release's Zenodo archive rather than in Git;
  - `checksum_only` if the terms allow neither. The table then stays out of Git (`.gitignore`) and anyone with the
    download can rebuild it.
- **Licence.** Check the licence before anything is built on it. Share-alike licences (CC BY-SA, CC BY-NC-SA, ODbL)
  limit what the release, which is CC BY-NC-SA 4.0, may carry from the source. Non-commercial terms are listed in the
  release README.
- **Reproducible.** Running the script again on the same downloads must give the same bytes.

## 2. Register it

Add a `Layer` to `LAYERS` in `scripts/release_layers.py`. The `CANCER` entry is the worked example. Every field you
fill is picked up by the build, the data dictionary, the workbook and the README; nothing else needs editing.

| Field | What it holds |
|---|---|
| `files`, `staging`, `licence`, `citation` | The staged files (`release_layers.py check` verifies them against MANIFEST at the start of every build), the scripts that make them, the terms, and how to cite it |
| `site_columns` and `site_values(sites, inputs)` | Site-table columns, as `(name, type, codes, description)`, and the function that fills them from the staged tables |
| `protein_columns` and `protein_values(proteins, inputs)` | The same for the protein table |
| `codes` | The code lists of coded columns (`data_dictionary.tsv` prints them) |
| `analyses` and `stats(sites, inputs, sequences, workdir)` | Statistics analyses with their description, and the function that returns their rows in the columns of `release_codes.STATISTICS_COLUMNS` |
| `supporting(sites, inputs, version)` and `supporting_rows` | Tables of the release's `supporting/` folder, as `{key: (path, frame, columns)}`, and what one row of each is (the README lists them) |
| `website(release, version, inputs, site_out)` | Writes the layer's website files from the release (the build's last step) |
| `topic` | Header colour of its site columns in the workbook |
| `folder` | Where its files are under `inputs/` (`layers` unless set) |
| `since`, `summary` and `changes(files)` | The release that adds it, what it adds, and sentences with its counts: `CHANGES.md` describes the layer under "New annotation layers" in that release |

A layer of knockouts, knockdowns or inhibitors needs no new code: add a `Perturbation` to `PERTURBATIONS` (its staged
peptide table, the columns of log2 change, p and adjusted p or `'bh'` (or `''` when the study gives no p), and the
column of the published call). Every arm is then called at twofold and q < 0.05 and by the published criteria, and
shown in the same tables and website cards. Set `primary='published'` (with `primary_why`) when the published call is
the one to show, as for 2 against 2 samples or a single experiment; `source_sites` when the study calls sites rather than
peptides; `localised` when some positions were assigned without the search localising them; `changed_only` when the
study lists its changed sites only.

Test the way the other tests are done, so that a layer's numbers mean what the rest of the release's do:

- **Sets of proteins** (classes, lists, locations): `enrichment.Job` with an explicit universe, which runs
  clusterProfiler `enricher`. Use `whole_universe=True` when the sets partition the background. The background is what
  the source measured (for example, the genes DepMap screened), not the whole proteome. See `cancer_stats`.
- **Residues** (a property of some arginines): compare methylated arginines with the other arginines of the same
  proteins.
  - Pooled: two-sided Fisher test, the odds ratio with 0.5 added to each cell.
  - Within protein: Mantel–Haenszel, each protein a stratum (`common.mantel_haenszel_ci`).
  - q is Benjamini–Hochberg within the analysis and method.
  - Sites the analyses leave out (column `counted` false) belong to neither group.
- **Whatever the source did not measure** is outside the test. A protein whose gene the source never measured is not
  a protein with a negative result.

## 3. Show it on the website

- **Per-protein data** go in 64 shard files named by `shard_of(accession)` (as `scripts/44c_export_depmap.py` and
  `scripts/36c` compute it). A page loads only its protein's file. A small layer (a few hundred kB) can write one file,
  which `loadLayer(name)` in `docs/assets/explorer.js` fetches once, when a card first needs it (`perturbations.json`,
  `reactome.json`, `paxdb.json`, `curation.json`).
- **The release's statistics of the layer** go into the layer's `meta.json` as `stats` rows with `stats_fields`. The
  website then shows the release's own numbers (`layerStats` in `docs/assets/explorer.js`) and does not compute them a
  second time.
- **Protein-set charts** are `hbars` specs. Rows hold `value` (log2 fold enrichment), `q`, `extra` and `title`. When a
  selection is active, `selSetTest` retests the chart against the same background, so the selection gets its own test
  and the whole dataset stays marked on each bar.
- **Other figures** go through
  `figBlock({ id, title, subText, figSvg, figTsv, figCss })`. That gives them the SVG, PNG and TSV exports every chart
  has. `figCss` maps the figure's CSS classes to the export colours.
- **Changes to `explorer.js` or `explorer.css`**: bump the `?v=` in `docs/index.html`.

## 4. Check it

1. `python scripts/release_layers.py check --inputs inputs`.
2. Build into scratch, not into `releases/` or `docs/data/`:

   ```
   RELEASE_OUT=/tmp/rel SITE_OUT=/tmp/site PREVIOUS=releases/vX.Y.Z bash hpc/build_release.sh vNEXT /tmp/work
   ```

3. Every column that was there before must match the previous release, and every new number must be checked
   independently. A short script on the staged tables is enough.
4. Serve the website with the new data (`python -m http.server` in a copy of `docs/` whose `data/` is `SITE_OUT`).
   Open the Overview with and without a selection, and a protein page. Try each export.
5. Merge through `develop` as `CONTRIBUTING.md` says. A layer that adds columns or statistics to the release is a new
   version.
