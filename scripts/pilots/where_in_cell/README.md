# Where in the cell: pilot (preview of v1.5.0, not part of a release)

Each map is tested against the proteins it places: methylated proteins (release v1.4.1) per compartment, one-sided
hypergeometric with fold enrichment, and the Mantel-Haenszel odds ratio within PaxDb abundance deciles of the same
proteins (`abundance_strata.py`). `build_preview.py` writes `docs/data/subcellular.json`, which the website's "Where in
the cell" card (preview) and the protein page's "Where the maps place it" read. The file carries each protein's PaxDb
decile among the proteins each map places, and the HPA locations of every HPA-profiled protein, so that the card repeats
the tests within the deciles for a selection (Mantel-Haenszel in the browser; with all methylated proteins as the
selection it gives the pilots' counts, odds ratios and confidence intervals).

| Map | File read | Source and licence |
|---|---|---|
| Human Cell Map (Go et al. 2021, HEK293) | `preys-v1.txt`, `not-significant-v1.txt` | https://humancellmap.org/resources/downloads (CC BY 4.0) |
| Organelle IP (Hein et al. 2025, HEK293T) | `data/labels/cluster_annotation_Dec6.csv` | github.com/czbiohub-sf/Organelle_IP_analyses_and_figures (BSD-3) |
| Chromatin probability (Kustatscher et al. 2014) | `embj0033-0648-sd8.xlsx` (Supplementary Table 1) | EMBO J 33:648, open access (CC BY) |
| Human Protein Atlas 25.1 | release statistics v1.4.1 (`hpa_location`) | already in the release (CC BY 4.0) |

The inputs are not in the repository yet; v1.5.0 would stage them into `inputs/` with checksums, as for every layer.

    python pilot_hcm.py REPO HCM_DOWNLOADS OUT
    python pilot_oip.py REPO ORGANELLE_IP_REPO OUT
    python pilot_chrom.py REPO embj0033-0648-sd8.xlsx OUT
    python build_preview.py REPO OUT docs/data/subcellular.json
