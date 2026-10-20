# Licence

Copyright © 2026 the PTM_Pipeline authors (Shechter Lab, Albert Einstein College of Medicine).

This repository holds two kinds of material under two licences:

| Material | Licence |
|---|---|
| Source code: `scripts/`, `hpc/`, `envs/`, `config/`, and the explorer's `docs/**/index.html`, `docs/**/assets/*.js` and `docs/**/assets/*.css` | [MIT](#code-mit-license) |
| Data, figures and text: the releases in `releases/`, the explorer's data (`docs/**/data/`) and written content, and `provenance/` | [CC BY-NC-SA 4.0](#data-and-content-cc-by-nc-sa-40) |
| The frozen inputs in `inputs/` | each its source's licence, listed in `inputs/MANIFEST.tsv` |

Third-party data keep their own licences (see [below](#third-party-data)).

## Code: MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Data and content: CC BY-NC-SA 4.0

The releases and the exported data, figures and text of the Methylarginine Site Explorer are licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). In short, you may share and adapt them for non-commercial purposes if you:
- give credit,
- indicate any changes, and
- release adaptations under the same licence.

Please cite DeAngelo JD, Silverstein AM, Shechter D. *Protein arginine methylation rewires molecular interactions.* Cell Chem Biol 33(10):1404–1421 (2026), [doi:10.1016/j.chembiol.2026.09.011](https://doi.org/10.1016/j.chembiol.2026.09.011). Also cite the original sources for any data you use.

### Why this licence

The site union and annotations include data distributed under share-alike, non-commercial terms:
- PhosphoSitePlus® (CC BY-NC-SA 3.0)
- DisGeNET (CC BY-NC-SA 4.0)
- the pepR-MS crosslink compendium (CC BY-NC)

CC BY-NC-SA 4.0 is the licence under which the combined data can be redistributed while respecting all of them. Creative Commons recommends against its licences for software, so the code is MIT.

## Third-party data

Data from other resources remain under their own terms. The explorer's **Documentation** page lists each source with its version, licence and citation. Those with non-commercial terms:

| Source | Licence |
|---|---|
| PhosphoSitePlus® (downloaded 22 May 2021, as compiled in Maron et al. 2021 Table S5) | CC BY-NC-SA 3.0. Acknowledge PhosphoSitePlus®, www.phosphosite.org, and cite Hornbeck et al. 2015, doi:10.1093/nar/gku1267. |
| DisGeNET | CC BY-NC-SA 4.0 |
| pepR-MS (Sha, Kuster & Trendel, bioRxiv 2025, doi:10.1101/2025.09.05.674461) | CC BY-NC |
| Trendel et al., Cell 188:4424 (2025), doi:10.1016/j.cell.2025.04.037 (DNA crosslinks and DNA interactome) | CC BY-NC 4.0 (open access article and supplement) |
| CD-CODE (Rostam et al., Nat Methods 20:673, 2023) | CC BY-NC 4.0 (stated on cd-code.org) |

Most other sources are CC0 or CC BY 4.0, including UniProt, STRING, Complex Portal, Gene Ontology, AlphaFold DB, AlphaMissense, ProtVar, the PDB, InterPro, Open Targets, ChIP-Atlas, IntAct, the Human Protein Atlas, RBDmap (Castello et al. 2016, CC BY 4.0) and the CARM1 knockout data (Shishkova et al. 2017, CC BY 4.0). Data reproduced from journal supplements are cited to their articles and remain subject to the publishers' terms. These include:
- Maron et al. 2021
- ProMetheusDB

CPLM 4.0 states "All Rights Reserved" and no licence. The current website (v1.1.0) does not use it; the published explorer kept in `docs/v1.0.0/` shows lysine positions derived from it. Its download and the tables made from it are not redistributed: `inputs/MANIFEST.tsv` records only their checksums.
