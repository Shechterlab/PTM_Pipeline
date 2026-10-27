# Acidic IDRs and the modifications that change charge: a proposal

A methylarginine is one side of an electrostatic contact. The other side is often an acidic stretch: the acidic tracts of
nucleophosmin that hold R-rich motifs in the nucleolus (Mitrea et al. 2016), the acidic IDRs of histone chaperones that
catch basic histones (Warren & Shechter 2017), the alternating charge blocks of DDX4, whose droplets arginine methylation
dissolves (Nott et al. 2015), and of MED1 (Lyons et al. 2023). The modifications that tune the acidic side are different
from those of the basic side: phosphorylation by the acidophilic kinase CK2 adds negative charge, TTLL glutamylation
lengthens acidic chains (NAP1, Regnard et al. 2000; Xenopus Npm2 and Nap1, Onikubo et al. 2015 and Lorton et al. 2024),
and ADP-ribosylation can sit on glutamate, aspartate or serine. On the basic side, acetylation removes a lysine's
charge and citrullination an arginine's. This proposal is for a parallel view of both sides, in two steps: a layer of
this repository first, a separate site later.

Proposed branch: `feature/acidic-idr`, from `develop` (`CONTRIBUTING.md`). A layer that adds columns and statistics is a
new version (v1.5.0).

## What there is now

**Here.** The other-PTM sites come from the iPTMnet-format master table (`inputs/sources/ptm_master_iptmnet_human.tsv`;
undated, no build script; `inputs/MANIFEST.tsv`): 119,603 phosphorylation rows, 5,124 lysine acylation, 698
ubiquitin/SUMO and 245 lysine methylation, no glutamylation, no ADP-ribosylation, and kinases only where the source
names an enzyme. The protein map draws charge (EMBOSS, 5-residue mean) and UniProt compositional bias, but nothing
calls an acidic tract and no test asks about one.

**In [Shechterlab/MicroDomains](https://github.com/Shechterlab/MicroDomains)** (private). The acidic pole of the dtIDR
score (`Analysis/microdomains/acidic_dtidr_notes.md`): domain-edge flanks scored for disorder, net negative charge and
glutamate in an acidic run. Histone chaperones (curated set, OR 3.0) and HMG proteins top it, the histone fold is
depleted, and CK2-attributed phosphorylation is enriched there only (OR 9.1; OR 5.0 among phosphorylated flanks), while
arginine methylation is depleted (OR 0.64). Three things found while reading it bear on this proposal:

1. The acidic score is still computed over the old, PTM-selected universe (`tables/acidic/acidic_dtidr_scores.tsv`:
   35,141 flanks on 3,522 proteins), not the full-proteome rebuild the basic score moved to (126,706 flanks,
   `tables/basic/all_domain_flanks.tsv`). ASF1A, ASF1B, DAXX, HIRA, NPM2, CHAF1A, NAP1L4, CGAS and TUBA1A have flanks
   in the full table and none in the acidic one, so caveat 1 of its notes (absent for lack of a domain annotation)
   describes the old universe, not the proteins.
2. Its PTM input, `Analysis/data/ptm_sites_in_domain_flanks.tsv`, is said to be written by
   `PTM_Pipeline/scripts/33_export_flank_ptm_annotation.py`. No commit, branch or tag of this repository has that file;
   `scripts/33_*` is the neurodegenerative disease figure. Its source and PTM-group columns are those of the iPTMnet
   table above.
3. The CK2 result rests on 42 attributed sites, and its §9 asks for a second attribution source and a glutamylation
   site list. Both can come from UniProt (below).

## A preview from the staged inputs

`scripts/52_preview_acidic_tracts.py` reads only the reference proteome, MobiDB disorder, PaxDb and the v1.4.1 sites. An
acidic tract is a run of at least 6 D/E (MicroDomains' prospective-list rule), counted when at least half of it is
disordered: 1,493 tracts on 1,149 of the 20,416 reviewed proteins.

| Test | Result (tracts of >= 6 D/E) |
|---|---|
| (A) Methylated proteins carry a disordered acidic tract | 467 of 4,439 (10.5%) against 682 of 15,977 other proteins (4.3%); within PaxDb abundance deciles, MH OR 2.46 (2.15-2.81) |
| (B) In cis: a methylarginine lies within 30 aa of a tract, among the disordered arginines of the same protein | 8.4% against 10.7%; within protein, MH OR 0.56 (0.42-0.74), rotation p 0.0005, 271 informative proteins |
| (B') The same without RG/GR arginines | MH OR 0.72 (0.50-1.03) within 30 aa; 0.55 (0.40-0.76) within 60 aa |

Across tract lengths of 5, 6 and 8 and distances of 10, 30, 60 and 100 aa, (A) is 2.33-2.46 and (B) 0.41-0.62 (rotation
p <= 0.03 every time); (B') is 0.48-0.80, and its interval crosses 1 in 4 of the 12 settings, all at 10-30 aa with
tracts of 6 or 8. So methylated proteins are acidic-tract proteins, and their methylarginines keep away from the acidic
tracts on the same chain: the two charges sit in separate blocks, part of which (not all) is the RG context of the
methylated blocks. This is the in-cis half of MicroDomains' last open question (§9.4, whether an acidic dtIDR on one
partner reads a basic microdomain on the other), which is in trans. Caveats: reported sites, with abundance controlled
but not study coverage; the tract rule is a convention, hence the sweep; disorder is MobiDB's consensus.

## Recommendation: a layer here first, a separate site later

Not a separate site yet:

- **The data are the bottleneck, not the viewer.** UniProt annotates polyglutamylation at 19 human sites (below), and a
  search found no proteome-scale table of glutamylation sites to stage (van Dijk et al. 2008 screened for substrate
  proteins). A site built now would mostly show sequence predictions.
- **The question needs both sides on one map.** The explorer already draws every reviewed protein with charge,
  composition, phosphorylation and acylation rows. An acidic-tract row and a row of glutamylation and ADP-ribosylation
  are a layer (`ADDING_A_SOURCE.md`), and (B) can only be seen with methylarginines and acidic tracts on one chain.
- **Staging belongs here.** MANIFEST, licences and checksums live in this repository, and MicroDomains reads an
  unversioned export of it. A released table fixes that.
- **MicroDomains is private and grant-facing.** A public segment-level site would publish its rankings before the lab
  means to.

A separate site makes sense when the unit changes from a residue to a segment: a dtIDR explorer whose rows are the
ranked flanks of both poles, linked to this explorer for the residues. Build it when the MicroDomains rankings are to be
public, as a second page in `docs/` (`docs/dtidr/`, sharing `explorer.js` and the data files) or as Pages of its own
repository.

## Work packages on `feature/acidic-idr`

1. **UniProt modified residues, all of them.** Stage every MOD_RES of the reference release (2026_03) the way
   `scripts/00_stage_uniprot_methylarginine.py` stages the methylarginines (CC BY 4.0, `git`): accession, position,
   residue, term, class (phospho, acetyl, glutamyl, glycyl, ADP-ribosyl, citrulline, ...), the kinase UniProt names
   ("Phosphoserine; by CK2"), evidence code and class. Counts in the current UniProt release (reviewed human, queried
   27 October 2026; 2026_03 to be staged):

   | Annotation | Sites | Proteins |
   |---|---|---|
   | Phosphorylation "by CK2" | 177 (137 with ECO:0000269) | 91 |
   | 5-glutamyl polyglutamate | 19 | 14 (12 tubulins, CGAS, KLF4) |
   | ADP-ribosyl glutamic acid | 38 | 17 |
   | ADP-ribosylserine | 73 | 39 |
   | Citrulline | 117 | 65 |
   | N6-acetyllysine | 4,355 | 1,710 |

   The CK2 rows are the second attribution source MicroDomains asks for. Keep the iPTMnet table beside it, as now.
2. **Glutamylation, curated.** UniProt has no non-tubulin chaperone. A table in `inputs/curation/` (with a checker, as
   `scripts/49`), kept apart from every statistic since curation is not a sample: NAP-1 E356/E357 and NAP-2 E347/E348
   (Regnard et al. 2000; both are glutamates of NAP1L1 and NAP1L4 in 2026_03, numbering to be checked against the paper),
   the substrates of van Dijk et al. 2008, CGAS (Xia et al. 2016; UniProt's CGAS E286 cites no paper), and the
   Xenopus Npm2 and Nap1 sites (Onikubo et al. 2015; Lorton et al. 2024) flagged as non-human, with the human positions
   only where an alignment places them.
3. **Acidic tracts, computed at build.** From the reference and MobiDB, no download: a supporting table of tracts
   (protein, start, end, length, disordered fraction, net charge), protein columns (number of tracts, longest run) and a
   site column (distance to the nearest disordered tract). The rule follows the decision below.
4. **Register a `Layer`** in `scripts/release_layers.py`: tests (A) and (B) with the rotation p, and, with the UniProt
   attribution, CK2 phosphorylation in acidic tracts against the other disordered S/T of the same proteins
   (within protein, as the release's residue tests are done).
5. **Website.** A map row of acidic tracts; one row for glutamylation, ADP-ribosylation and citrullination; CK2 marked
   on the phosphorylation row; the layer's statistics in its `meta.json`. Bump `?v=` in `docs/index.html`.
6. **For MicroDomains.** The other-PTM sites, with source, enzyme and evidence, as a supporting table of the release, so
   MicroDomains reads a version rather than an export; then rescore its acidic pole on the full-proteome flanks.

## Decisions for the lab

- The tract rule: a run of >= 6 D/E, a net-charge window, or MicroDomains' acidic score on flanks.
- Whether non-human glutamylation sites enter the curation table, and how they are placed.
- Whether UniProt replaces the iPTMnet table for crosstalk on the website or sits beside it.
- When the dtIDR site goes public.

## References

- Lorton et al., iScience 27:109458 (2024): [Einstein research profile](https://einstein.elsevierpure.com/en/publications/glutamylation-of-npm2-and-nap1-acidic-disordered-regions-increase)
- Lyons et al., Cell 186:327 (2023), [doi:10.1016/j.cell.2022.12.013](https://doi.org/10.1016/j.cell.2022.12.013)
- Mitrea et al., eLife 5:e13571 (2016), [doi:10.7554/eLife.13571](https://doi.org/10.7554/eLife.13571)
- Nott et al., Mol Cell (2015), [PMC4352761](https://pmc.ncbi.nlm.nih.gov/articles/PMC4352761/)
- Onikubo et al., Cell Rep (2015), [PMC4567554](https://pmc.ncbi.nlm.nih.gov/articles/PMC4567554/)
- Regnard et al., J Biol Chem 275:15969 (2000), [doi:10.1074/jbc.M000045200](https://doi.org/10.1074/jbc.M000045200)
- van Dijk et al., J Biol Chem 283:3915 (2008); reviewed in Janke, Rogowski & van Dijk, EMBO Rep 9:636 (2008),
  [PMC2475320](https://pmc.ncbi.nlm.nih.gov/articles/PMC2475320/)
- Warren & Shechter, J Mol Biol 429:2401 (2017), [doi:10.1016/j.jmb.2017.06.005](https://doi.org/10.1016/j.jmb.2017.06.005)
- Xia et al., Nat Immunol (2016), [doi:10.1038/ni.3356](https://doi.org/10.1038/ni.3356)
