# Studies and databases behind the sites, v1.2.0

Every site rests on records: rows of the source tables, each crediting the study or database that reported the site. A record is placed on a reviewed UniProt entry by its own sequence evidence (a window centred on the site, or a peptide with or without the site marked), and only by its reported position when it gives no sequence; every placed record is then checked against that evidence, and every site is an arginine of the reference sequence. The records of a study read from its own table add a site only when the study localised the methyl group to one arginine and the peptide (or, without one, the window) names exactly one arginine in the whole reviewed proteome, in a protein the record names; the rows of the compiled tables (Maron et al. 2021 Table S5, the PTM master table) carry no localisation and add none.

Each section gives where the study's records come from, what they hold, what became of them, and records traced from the source table to the reviewed proteome: a placed record, with its evidence aligned to the reference sequence (brackets mark the site), and a record left out, with the reason. The placement codes are explained in `data_dictionary.tsv`; every record is in `supporting/source_records`.

| Study or database | Records from | Records | Placed | Sites | Reported by no other source |
|---|---|---:|---:|---:|---:|
| [Larsen 2016](#larsen-2016) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Larsen et al. 2016 Table S1, staged (scripts/00_stage_study_sites.py) | 17,766 | 7,297 | 7,260 | 1,645 |
| [PhosphoSitePlus (PSP)](#phosphositeplus-psp) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies | 8,375 | 7,900 | 7,443 | 1,232 |
| [Maron 2021](#maron-2021) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Maron et al. 2021 Tables S2 and S3 with the authors' Proteome Discoverer and Byonic search results, staged (scripts/00_stage_study_sites.py) | 6,108 | 3,749 | 1,482 | 372 |
| [Guo 2014](#guo-2014) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Guo et al. 2014 supplementary Tables S2 and S3, staged (scripts/00_stage_study_sites.py) | 4,729 | 2,999 | 1,569 | 48 |
| [Musiani 2019](#musiani-2019) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Musiani et al. 2019 Tables S1 and S2, staged (scripts/00_stage_study_sites.py) | 4,223 | 2,093 | 1,299 | 226 |
| [Fong 2019](#fong-2019) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Fong et al. 2019 supplementary Table S4, staged (scripts/00_stage_study_sites.py) | 3,841 | 2,185 | 1,270 | 173 |
| [ProMetheusDB](#prometheusdb) | ProMetheusDB (Massignani et al. 2022) supplementary table mmc2 | 1,598 | 292 | 281 | 77 |
| [Li 2021](#li-2021) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Li et al. 2021 Supplementary Data 1, 3, 5 and 6, staged (scripts/00_stage_study_sites.py) | 2,880 | 329 | 319 | 43 |
| [Fedoriw 2019](#fedoriw-2019) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Fedoriw et al. 2019 supplementary Tables S5 and S7, staged (scripts/00_stage_study_sites.py) | 3,610 | 536 | 409 | 31 |
| [UniProt](#uniprot) | the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records); UniProt release 2026_03 MOD_RES with experimental evidence, staged (scripts/00_stage_study_sites.py) | 1,016 | 254 | 212 | 44 |
| [Lim 2020](#lim-2020) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Lim et al. 2020 Additional files 1 and 2, staged (scripts/00_stage_study_sites.py) | 2,001 | 1,254 | 561 | 86 |
| [Wei 2020](#wei-2020) | Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Wei et al. 2020 Table 1, staged (scripts/00_stage_study_sites.py) | 578 | 283 | 283 | 59 |
| [IEDB](#iedb) | the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records); IEDB epitopes with a methylated arginine (query API), staged (scripts/00_stage_study_sites.py) | 166 | 43 | 43 | 20 |
| [Hartel 2019](#hartel-2019) | Hartel et al. 2019 supplementary Tables S1-S15, staged (scripts/00_stage_study_sites.py) | 2,806 | 2,684 | 1,447 | 125 |
| [Geoghegan 2015](#geoghegan-2015) | Geoghegan et al. 2015 Supplementary Table 1, staged (scripts/00_stage_study_sites.py) | 4,536 | 4,069 | 2,233 | 203 |
| [Sylvestersen 2014](#sylvestersen-2014) | Sylvestersen et al. 2014 supplementary tables, staged (scripts/00_stage_study_sites.py) | 1,281 | 993 | 889 | 29 |
| [Lu 2025](#lu-2025) | Lu et al. 2025 Dataset EV1, staged (scripts/00_stage_study_sites.py) | 743 | 155 | 155 | 11 |
| [Ma 2023](#ma-2023) | Ma et al. 2023 Table S1 (PRIDE PXD038203), staged (scripts/00_stage_study_sites.py) | 403 | 403 | 403 | 46 |

## Larsen 2016

Larsen SC et al., Sci Signal 9:rs9 (2016), [doi:10.1126/scisignal.aaf7329](https://doi.org/10.1126/scisignal.aaf7329).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Larsen et al. 2016 Table S1, staged (scripts/00_stage_study_sites.py). **Terms:** AAAS (Sci Signal supplementary Table S1: derived records only, the table is not redistributed).

> HEK293 cells, peptides fractionated at high pH and enriched with methylarginine antibodies, MaxQuant with peptide, protein and site FDR 0.01, monomethylation of R and K searched. Read from the article's Table S1 (8,900 sites, each with the leading protein, the position, the localisation probability and a 31-residue window): arginine sites at localisation probability >= 0.75. The leading protein is often an isoform or an unreviewed entry, so the site is placed from its window.

17,766 records: 7,297 placed on 7,260 sites (7,260 of them in the previous version, 1,645 reported by no other source). Not placed: 8,884 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 1,244 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 281 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 54 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)); 6 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome).

```text
record 28965: reported A0AV96 at R332 (localisation probability 1, 1 source row, tables Table S1)
placed on A0AV96:R332 (sequence_unique_canonical; evidence check: agrees); the window of 31 residues occurs once in the reviewed proteome
evidence              KPVDKEQYSRYQKAA[R]GGGAAEAAQQPSYVY
reference   311 LEVTLAKPVDKEQYSRYQKAA[R]GGGAAEAAQQPSYVYSCDPYT 353
```

```text
record 29685: reported G3V2F7 at R348, window of 31 residues KWQNSYSIKVVLQELRRLMMSKENMKLPQPP
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 2 sites in the reviewed proteome: Q13404:R125, Q15819:R123
```

```text
record 28984: reported A3KN83 at R1024, window of 31 residues VAKRLESLGALTHGDRRATESRDLSKYNFEN
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: Q9Y2G9:R884
```

## PhosphoSitePlus (PSP)

Hornbeck PV et al., Nucleic Acids Res 43:D512-D520 (2015), [doi:10.1093/nar/gku1267](https://doi.org/10.1093/nar/gku1267). Downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5.

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies. **Terms:** CC BY-NC-SA 3.0 Unported (PhosphoSitePlus terms of 2021, inputs/licences/), compiled in Maron 2021 Table S5 (CC BY-NC-ND 4.0, open-access supplement).

> Curated database, trusted as curated. Through Maron Table S5: the sites PhosphoSitePlus listed on 22 May 2021, each with a 15-residue window and position, without the study behind it, a localisation score or a degree. A record is placed where its window names exactly one arginine in the reviewed proteome, in an entry the record names.

8,375 records: 7,900 placed on 7,443 sites (7,443 of them in the previous version, 1,232 reported by no other source). Not placed: 436 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 26 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 13 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)).

```text
record 13126: reported A0AV96 at R332
placed on A0AV96:R332 (sequence_unique_canonical; evidence check: agrees); the window of 15 residues occurs once in the reviewed proteome
evidence              SRYQKAA[R]GGGAAEA
reference   319 VDKEQYSRYQKAA[R]GGGAAEAAQQPSY 345
```

```text
record 13138: reported A0FGR8 at R23, window of 15 residues VPSEGRWRGAEPPGI
not placed (unresolved_study_not_found): a localised record whose evidence names no arginine anywhere in the reference proteome
the evidence occurs nowhere in the reviewed proteome
```

## Maron 2021

Maron MI et al., iScience 24:102971 (2021), [doi:10.1016/j.isci.2021.102971](https://doi.org/10.1016/j.isci.2021.102971). The study's own data: Tables S2 and S3, localised with the authors' Proteome Discoverer and Byonic search results.

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Maron et al. 2021 Tables S2 and S3 with the authors' Proteome Discoverer and Byonic search results, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY-NC-ND 4.0 (open-access supplement, Tables S2 and S3), with localisation from the authors' search results (derived records only, the files are not redistributed).

> A549 cells treated with GSK591 or MS023, sequential Rme1, Rme2a and Rme2s antibody immunoprecipitation, Byonic in Proteome Discoverer 2.4 with mono- and dimethylation of R and K searched, 1 % FDR, trypsin and GluC digests. Tables S2 and S3 give no localisation score, the authors' search results of the same immunoprecipitations do: Proteome Discoverer peptide groups with ptmRS site probabilities (a residue named only at >= 75 %), and the Byonic PSMs of the GluC and trypsin searches with their Delta Mod score. A mark of the tables is kept where the search localises a methyl group to that arginine of the same peptide (ptmRS >= 75 %, or Delta Mod > 20 in every high-confidence PSM of that peptide and set of marks above it, with the search's degree when it differs, noted), or where its peptide holds no other arginine or lysine. Marks the confident PSMs disagree on are kept as lower-confidence evidence.

6,108 records: 3,749 placed on 1,482 sites (1,482 of them in the previous version, 372 reported by no other source). Not placed: 1,829 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 506 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 18 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 6 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome).

```text
record 42073: reported Q9Y4Z0 at R125 (289 source rows, tables Table S2, Table S3, GluC peptide, Table S2 gives me1, the search localises me2, Byonic Delta Mod 108.9, Table S3 gives me1, the search localises me2)
placed on Q9Y4Z0:R125 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             VVAKGRGRGGLQQQKQQKGRGMGGAGRGVFGGRGRGGIPGTG[R]GQPE
reference   77 DMVKEEVVAKGRGRGGLQQQKQQKGRGMGGAGRGVFGGRGRGGIPGTG[R]GQPEKKPGRQ 135
```

```text
record 37882: reported C9JG80 at R1048, marked peptide R(di)RVGDVEPSR
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 17 sites in the reviewed proteome: A6NJU9:R1048, A6NJU9:R1077, A6NJU9:R1106, A8MRT5:R1043, A8MRT5:R1072, A8MRT5:R1101, ...
```

## Guo 2014

Guo A et al., Mol Cell Proteomics 13:372-387 (2014), [doi:10.1074/mcp.O113.027870](https://doi.org/10.1074/mcp.O113.027870).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Guo et al. 2014 supplementary Tables S2 and S3, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (supplementary Tables S2 and S3).

> HCT116 cells, monomethylarginine (Me-R4-100, D5A12) and ADMA antibodies, SEQUEST with a 1 % PSM FDR, mono- and dimethylation of R and K searched. Read from Tables S2 (peptide-spectrum matches with their ModScore) and S3 (sites with a 13-residue window): a site is kept where a PSM matching its window localises it at ModScore > 13, the paper's threshold, the degree read from the precursor mass. The positions number RefSeq proteins, so the record is placed from its peptide.

4,729 records: 2,999 placed on 1,569 sites (1,566 of them in the previous version, 48 reported by no other source). Not placed: 1,471 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 169 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 58 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 27 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)); 5 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome).

```text
record 42468: reported P05787 at R47 (111 source rows, tables Table S2, Table S3, ModScore 453.585, RefSeq NP_002264)
placed on P05787:R47 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             VGSSNF[R]GGLGGGYGGASGMGGITAVTVNQSLLSPLVLEVDPNIQAVR
reference   35 SSSFSRVGSSNF[R]GGLGGGYGGASGMGGITAVTVNQSLLSPLVLEVDPNIQAVRTQEKEQ 94
```

```text
record 44466: reported Q92804 at R456, marked peptide GGGYGGDR(me)GGGYGGD
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 4 sites in the reviewed proteome: Q92804:R459, Q92804:R467, Q92804:R475, Q92804:R562
```

```text
record 42717: reported P19338 at R660, marked peptide GR(me)GGFGGRGGFR
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: P19338:R681
```

## Musiani 2019

Musiani D et al., Sci Signal 12:eaat8388 (2019), [doi:10.1126/scisignal.aat8388](https://doi.org/10.1126/scisignal.aat8388).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Musiani et al. 2019 Tables S1 and S2, staged (scripts/00_stage_study_sites.py). **Terms:** AAAS (Sci Signal supplementary Tables S1 and S2: derived records only, the tables are not redistributed).

> HeLa cells, SILAC with GSK591, SDMA then monomethylarginine antibody enrichment and the input, MaxQuant (1.6.1.0 in the PRIDE PXD009070 deposit) with a 1 % FDR, and heavy-methyl SILAC for the orthogonal validation. Read from the article's Table S1 (every site and identification, with the leading protein, position, 15-residue window and localisation probability) and Table S2 (the heavy-methyl SILAC pairs, Spreadsheet 2 naming the gene only, so its entry is the gene's reviewed entry that holds the peptide once): arginine sites at localisation probability >= 0.75, the best of the heavy and light peptides in Table S2. Where the window centres on more than one marked arginine (repeats such as the TAF15 RGG repeats), the record is placed from its window alone, and a window repeated in the protein adds no site.

4,223 records: 2,093 placed on 1,299 sites (1,283 of them in the previous version, 226 reported by no other source). Not placed: 1,712 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 324 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 56 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 31 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 7 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)).

```text
record 46040: reported P09651; H0YH80 at R206 (localisation probability 1, 31 source rows, tables Table S1, Table S2 Spreadsheet 1, Table S2 Spreadsheet 2, SILAC Rev_Sym, SILAC Rev_Input, SILAC Fwd_Input, SILAC Fwd_Sym, SILAC Rev_Mono, SILAC Fwd_Mono, LcPrb > 0.75 (as given))
placed on P09651:R206 (sequence_paralog_named; evidence check: one_of_several); the marked peptide occurs 2 times in the reviewed proteome
evidence              SGSGNFGGG[R]GGGFGGNDNFGRGGNFSGRGGFGGSR
reference   191 SSQRGRSGSGNFGGG[R]GGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYG 238
```

```text
record 45464: reported A0A1C7CYW6 at R24, marked peptide VIVGTIR(me)PR
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 2 sites in the reviewed proteome: P0DI83:R24, P0DI83:R33
```

```text
record 45486: reported A0A1W2PQP7 at R59, marked peptide TLR(di)GWSRSSR
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: Q6UXB0:R59
```

## Fong 2019

Fong JY et al., Cancer Cell 36:194-209 (2019), [doi:10.1016/j.ccell.2019.07.003](https://doi.org/10.1016/j.ccell.2019.07.003).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Fong et al. 2019 supplementary Table S4, staged (scripts/00_stage_study_sites.py). **Terms:** Elsevier user licence (supplementary Table S4: derived records only, the table is not redistributed), PRIDE PXD012007 CC0.

> NB4 acute promyelocytic leukaemia cells, SDMA and monomethylarginine (and, with MS023, ADMA) antibody enrichment, MaxQuant with a 1 % FDR. Read from Table S4: arginine sites at localisation probability >= 0.75 and Andromeda score >= 25 (the paper accepted localisation >= 0.5), with the marked peptide and degree.

3,841 records: 2,185 placed on 1,270 sites (1,256 of them in the previous version, 173 reported by no other source). Not placed: 1,298 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 322 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 19 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 16 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 1 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)).

```text
record 48746: reported P38159 at R125 (localisation probability 1, 348 source rows, tables Table S4)
placed on P38159:R125 (sequence_paralog_named; evidence check: one_of_several); the marked peptide occurs 2 times in the reviewed proteome
evidence              GPPS[R]GGHMDDGGYSMNFNMSSSR
reference   115 GSGGTRGPPS[R]GGHMDDGGYSMNFNMSSSRGPLPVK 150
```

```text
record 48544: reported P18583 at R1039, marked peptide SMMSAYER(me)
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 4 sites in the reviewed proteome: P18583:R1039, P18583:R1047, P18583:R1064, P18583:R1072
```

```text
record 47922: reported A0A0A6YYC7 at R1, marked peptide GR(me)AAAAAAAAAVSR
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: Q96JP5:R58
```

## ProMetheusDB

Massignani E et al., Mol Cell Proteomics 21:100243 (2022), [doi:10.1016/j.mcpro.2022.100243](https://doi.org/10.1016/j.mcpro.2022.100243).

**Records from** ProMetheusDB (Massignani et al. 2022) supplementary table mmc2. **Terms:** CC BY-NC-ND 4.0 (open-access supplement mmc2).

> Methyl-proteome table mmc2 of the heavy-methyl SILAC and SILAC experiments behind ProMetheusDB (SK-OV-3, HeLa, NB4 and U2OS cells), MaxQuant with mono- and dimethylation of R and K and trimethylation of K searched: marked peptides, 15-residue window and degree (me or di) of each site, no localisation score. A record is kept where one of its peptides holds no other arginine or lysine, or holds no lysine and has every arginine dimethylated, and is placed from those peptides.

1,598 records: 292 placed on 281 sites (281 of them in the previous version, 77 reported by no other source). Not placed: 1,236 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 69 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 1 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome).

```text
record 27365: reported P24666 at R28
placed on P24666:R28 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             SPIAEAVF[R]
reference   14 LGNICRSPIAEAVF[R]KLVTDQ 34
```

```text
record 27363: reported Q9UKV3 at R14, marked peptide TSGGTR(di)GVLSGNR
not placed (unresolved_not_localised): the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)
the evidence names 1 site in the reviewed proteome: Q9UKV3:R14
```

## Li 2021

Li WJ et al., Nat Commun 12:1946 (2021), [doi:10.1038/s41467-021-21963-1](https://doi.org/10.1038/s41467-021-21963-1).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Li et al. 2021 Supplementary Data 1, 3, 5 and 6, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (Supplementary Data 1, 3, 5 and 6).

> HEK293 cells (SILAC) after PRMT4, PRMT5 or PRMT7 knockdown or SGC3027 treatment, monomethylarginine then ADMA or SDMA antibody enrichment, Proteome Discoverer 2.2 with FDR < 1 %, methylation of R and K searched, no localisation score given. Read from Supplementary Data 1, 3, 5 and 6, the sites whose methylation fell at least two-fold or vanished: a site is kept where its peptide holds no other arginine or lysine.

2,880 records: 329 placed on 319 sites (317 of them in the previous version, 43 reported by no other source). Not placed: 1,486 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 1,059 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 5 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 1 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome).

```text
record 50634: reported P04264 at R82 (6 source rows, tables PRMT4-regulated me1 (FC>2), PRMT4-regulated me1 (abolished), PRMT5-regulated me1 (FC>2), PRMT5-regulated me1 (abolished), PRMT7-regulated me1 (FC>2), PRMT7-regulated me1 (abolished), PRMT7 knockdown, PRMT4 knockdown, PRMT5 knockdown)
placed on P04264:R82 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             SISISVA[R]
reference   69 NLGGSKSISISVA[R]GGGRGS 88
```

```text
record 50811: reported P18583 at R1031, marked peptide SMMSPMAER(me)
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 2 sites in the reviewed proteome: P18583:R1031, P18583:R1056
```

## Fedoriw 2019

Fedoriw A et al., Cancer Cell 36:100-114 (2019), [doi:10.1016/j.ccell.2019.05.014](https://doi.org/10.1016/j.ccell.2019.05.014).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Fedoriw et al. 2019 supplementary Tables S5 and S7, staged (scripts/00_stage_study_sites.py). **Terms:** Elsevier user licence (supplementary Tables S5 and S7: derived records only, the tables are not redistributed).

> DLBCL and pancreatic cancer cell lines with a type I PRMT inhibitor, monomethylarginine, ADMA and SDMA antibody enrichment, SEQUEST with a 2.5 % FDR, no localisation score given. Read from Tables S5 and S7 (gene symbol and marked peptide, the peptides whose methylation changed): a site is kept where its peptide holds one arginine, or every arginine dimethylated, and occurs once in the entries of its gene, which gives the position.

3,610 records: 536 placed on 409 sites (399 of them in the previous version, 31 reported by no other source). Not placed: 2,105 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 951 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 18 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome).

```text
record 54367: reported Q96AE4 at R359 (9 source rows, tables Table S5, SDMA immunoprecipitation)
placed on Q96AE4:R359 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence              SVQAGNPGGPGPGG[R]GRGRGQGNWNMGPPGGLQEFNFIVPTGK
reference   339 ITDLLRSVQAGNPGGPGPGG[R]GRGRGQGNWNMGPPGGLQEFNFIVPTGKTGLIIG 393
```

## UniProt

UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025), [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010).

**Records from** the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records); UniProt release 2026_03 MOD_RES with experimental evidence, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (UniProt release 2026_03 MOD_RES).

> UniProt release 2026_03 annotations of methylated arginines (MOD_RES) with experimental evidence from the literature (ECO:0000269), curated residue by residue. Annotations from large-scale proteomics (ECO:0007744), which repeat the studies read here, and those by similarity are left out. Position only, placed where the residue is an arginine of the entry.

1,016 records: 254 placed on 212 sites (201 of them in the previous version, 44 reported by no other source). Not placed: 762 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation).

```text
record 54923: reported O00267 at R698 (2 source rows, tables MOD_RES, asymmetric (UniProt), symmetric (UniProt))
placed on O00267:R698 (same_accession_position_verified): no sequence, the reported residue is an arginine
reference   691 SGGMSRG[R]GRRDNEL 705
```

## Lim 2020

Lim Y et al., Proteome Sci 18:6 (2020), [doi:10.1186/s12953-020-00162-8](https://doi.org/10.1186/s12953-020-00162-8).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Lim et al. 2020 Additional files 1 and 2, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (Additional files 1 and 2).

> Pooled colorectal cancer tissue from 10 patients, monomethylarginine and ADMA antibody enrichment, SEQUEST with a 2.5 % FDR. Read from Tables S1 and S2: sites with ModScore > 13 (the paper's threshold), each with peptide, window and degree. The tables write some leucines as I: L is written back where the reported entry, or one site of the reviewed proteome, matches with I and L taken as one residue.

2,001 records: 1,254 placed on 561 sites (560 of them in the previous version, 86 reported by no other source). Not placed: 651 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 59 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 29 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 5 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)); 3 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome).

```text
record 55439: reported P05787 at R47 (119 source rows, tables Table S1, ModScore 1000)
placed on P05787:R47 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             VGSSNF[R]GGLGGGYGGASGMGGITAV
reference   35 SSSFSRVGSSNF[R]GGLGGGYGGASGMGGITAVTVNQSL 72
```

```text
record 56225: reported Q92804 at R459, marked peptide GGGYGGDR(me)GGGYGGDR
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 4 sites in the reviewed proteome: Q92804:R459, Q92804:R467, Q92804:R475, Q92804:R562
```

```text
record 55535: reported P09651 at R225, marked peptide GGNFSGR(me)GGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMK
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: A0A2R8Y4L2:R225
```

## Wei 2020

Wei M et al., Front Cell Dev Biol 8:545934 (2020), [doi:10.3389/fcell.2020.545934](https://doi.org/10.3389/fcell.2020.545934).

**Records from** Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies; Wei et al. 2020 Table 1, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (Table 1).

> HPDE6c7 and PANC-1 cells, ADMA antibody enrichment, MaxQuant 1.6.0.16 with FDR < 0.01, Andromeda score > 40 and localisation probability > 0.75. Read from Table 1: the 289 sites with peptide and 31-residue window, 137 of them numbered on an isoform, placed from the peptide.

578 records: 283 placed on 283 sites (281 of them in the previous version, 59 reported by no other source). Not placed: 289 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 6 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome).

```text
record 56522: reported A6NGB9 at R46 (localisation probability 1, 1 source row, tables Table 1)
placed on A6NGB9:R46 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence             G[R]SALLADIQQGTR
reference   39 RRADPKG[R]SALLADIQQGTRLRKVTQ 64
```

## IEDB

Vita R et al., Nucleic Acids Res 47:D339-D343 (2019), [doi:10.1093/nar/gky1006](https://doi.org/10.1093/nar/gky1006).

**Records from** the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records); IEDB epitopes with a methylated arginine (query API), staged (scripts/00_stage_study_sites.py). **Terms:** IEDB query API, terms at www.iedb.org/terms_of_use_v3.php (epitopes the iPTMnet-format master table cites).

> Immune Epitope Database epitopes (HLA-presented peptides) with a methylated arginine, the 85 the iPTMnet-format master table cites, read from the IEDB query API: a site is kept where the epitope holds no other arginine or lysine. The degree is not given.

166 records: 43 placed on 43 sites (43 of them in the previous version, 20 reported by no other source). Not placed: 85 `superseded_by_study_table` (a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation); 38 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)).

```text
record 56812: reported A7MCY6 at R427 (1 source row, tables IEDB epitope, IEDB epitope 1324010)
placed on A7MCY6:R427 (sequence_paralog_named; evidence check: one_of_several); the marked peptide occurs 11 times in the reviewed proteome
evidence              [R]PPPPPPP
reference   421 CPAPQP[R]PPPPPPPGERTLA 440
```

## Hartel 2019

Hartel NG et al., Mol Cell Proteomics 18:2149-2164 (2019), [doi:10.1074/mcp.RA119.001625](https://doi.org/10.1074/mcp.RA119.001625). HEK293T and LN229 cells. ADMA and SDMA called from characteristic neutral losses (dimethylamine, methylamine) of dimethylarginine peptides (reported_methyl_state Rme2a, Rme2s), the other dimethylarginines as Rme2.

**Records from** Hartel et al. 2019 supplementary Tables S1-S15, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (the article and its supplement, from 2020-09-15, and PRIDE PXD012357, CC0), records staged from supplementary Tables S1-S15.

> Peptide-spectrum matches and peptide groups (Proteome Discoverer) with the methylated arginine marked, Percolator q <= 0.01 and ptmRS localisation >= 75 %. LN229 Tables S1-S2 give no positions: their rows take the position Tables S3-S15 give for the same protein, peptide and residue, and are dropped without one.

2,806 records: 2,684 placed on 1,447 sites (1,339 of them in the previous version, 125 reported by no other source). Not placed: 89 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 19 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)); 14 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome).

```text
record 58064: reported Q00839 at R739 (neutral loss ADMA, localisation probability 1, 864 source rows, tables S1, S2, S3, S4, S5, S6, S9, S11, S12, S9 PRMT1 KD log2 -0.44 ns, S9 PRMT1 KD log2 2.77 ns, S11 PRMT1 KD log2 -0.21 ns, S12 PRMT1 KD log2 -0.34 ns)
placed on Q00839:R739 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence              RGNMPQ[R]GGGGGGSGGIGYPYPR
reference   727 RGGYNRRGNMPQ[R]GGGGGGSGGIGYPYPRAPVFPG 761
```

```text
record 59032: reported Q92804 at R415, marked peptide GGDR(me)GGYGGDR
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 6 sites in the reviewed proteome: Q92804:R415, Q92804:R483, Q92804:R498, Q92804:R505, Q92804:R512, Q92804:R519
```

```text
record 59076: reported Q92843-2 at R278, marked peptide TNR(di)PGISTTDRGFPRAR
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: Q86U42:R251
```

## Geoghegan 2015

Geoghegan V et al., Nat Commun 6:6758 (2015), [doi:10.1038/ncomms7758](https://doi.org/10.1038/ncomms7758). Jurkat and primary human T cells, heavy-methyl SILAC with methylarginine antibodies.

**Records from** Geoghegan et al. 2015 Supplementary Table 1, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (the article and its Supplementary Table 1).

> Heavy-methyl SILAC peptide pairs with the methylated arginine marked, q <= 0.05 (0.93 % FDR for the paired set), localised at ModLS > 8 or unambiguously (rows at ModLS 8 or below are kept as lower-confidence evidence). Positions number the first accession of the protein group, often an unreviewed entry, so the site is placed from the peptide.

4,536 records: 4,069 placed on 2,233 sites (2,065 of them in the previous version, 203 reported by no other source). Not placed: 235 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 181 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 43 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 8 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)).

```text
record 60051: reported B4DHY1; P31942; P31942-2; P31942-3; P31942-4 at R235 (59 source rows, tables ST1, ModLS unambiguous)
placed on P31942:R343 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence              MGMGNNYSGGYGTPDGLGGYGRGGGGSGGYYGQGGMSGGGW[R]GMY
reference   296 YGSVGRMGMGNNYSGGYGTPDGLGGYGRGGGGSGGYYGQGGMSGGGW[R]GMY 346
```

```text
record 61903: reported K7EKB2; Q92804; K7EN07; Q92804-2; K7EPT6 at R80, marked peptide GGDR(me)GGY
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 7 sites in the reviewed proteome: Q92804:R415, Q92804:R483, Q92804:R498, Q92804:R505, Q92804:R512, Q92804:R519, ...
```

```text
record 61210: reported F8W646; H0YH80; P09651-2 at R116, marked peptide GGSR(me)GGGGYGGSGDGYNGFGNDGSNF
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: A0A2R8Y4L2:R232
```

## Sylvestersen 2014

Sylvestersen KB et al., Mol Cell Proteomics 13:2072-2088 (2014), [doi:10.1074/mcp.O113.032748](https://doi.org/10.1074/mcp.O113.032748). HEK293T cells, monomethylarginine antibody enrichment.

**Records from** Sylvestersen et al. 2014 supplementary tables, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (the article and its supplement, from 2020-10-12, and PRIDE PXD000559, CC0).

> MaxQuant site tables with a 13-residue window and the peptide, localisation probability >= 0.75. The export writes every leucine as isoleucine: L is written back where the reported entry, or else one site of the reviewed proteome, matches with I and L taken as one.

1,281 records: 993 placed on 889 sites (869 of them in the previous version, 29 reported by no other source). Not placed: 175 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)); 96 `unresolved_study_not_found` (a localised record whose evidence names no arginine anywhere in the reference proteome); 10 `unresolved_study_not_unique` (a localised record whose evidence names more than one site in the reference proteome); 7 `unresolved_study_other_entry` (a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)).

```text
record 64235: reported A0A5E8 at R642 (localisation probability 0.9446, 2 source rows, tables 032748-4 All identified MMA sites, 032748-4 Triplicate experiment - MMA)
placed on Q99501:R642 (sequence_unique_candidate; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence              MDTQPD[R]KPSR
reference   630 TWARGRMDTQPD[R]KPSRIPTPRG 652
```

```text
record 64811: reported Q07283 at R1404, marked peptide QER(me)EQQLR
not placed (unresolved_study_not_unique): a localised record whose evidence names more than one site in the reference proteome
the evidence names 2 sites in the reviewed proteome: Q07283:R1404, Q07283:R1657
```

```text
record 64253: reported A4PIV7 at R9, marked peptide SVAFAAPR(me)
not placed (unresolved_study_other_entry): a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)
the evidence names 1 site in the reviewed proteome: Q15532:R9
```

## Lu 2025

Lu L et al., EMBO Rep 26:5649-5672 (2025), [doi:10.1038/s44319-025-00599-y](https://doi.org/10.1038/s44319-025-00599-y). HeLa cells, enrichment with the SMN Tudor domain.

**Records from** Lu et al. 2025 Dataset EV1, staged (scripts/00_stage_study_sites.py). **Terms:** CC BY 4.0 (the article and Dataset EV1).

> The 743 dimethylarginine sites of Dataset EV1 come without peptides or localisation scores (SEQUEST HT, Percolator 1 % FDR). A site is kept only where a PSM localises it by mass, every arginine of its peptide dimethylated (155 sites), and the other 588 are kept as lower-confidence evidence.

743 records: 155 placed on 155 sites (145 of them in the previous version, 11 reported by no other source). Not placed: 588 `unresolved_not_localised` (the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)).

```text
record 65535: reported O14979 at R408 (1 source row, tables EV1 Arginine dimethylation sites, EV1 PSMs, 5 PSMs with every arginine dimethylated)
placed on O14979:R408 (sequence_unique_canonical; evidence check: agrees); the marked peptide occurs once in the reviewed proteome
evidence              AS[R]GGGNHQNNYQPY
reference   400 QSTYGKAS[R]GGGNHQNNYQPY 420
```

## Ma 2023

Ma M et al., J Am Soc Mass Spectrom 34:1692-1700 (2023), [doi:10.1021/jasms.3c00154](https://doi.org/10.1021/jasms.3c00154). Breast cancer patient-derived xenografts grown in mice, ADMA antibody enrichment. A peptide assigned to a human protein can come from the mouse stroma where the two sequences agree.

**Records from** Ma et al. 2023 Table S1 (PRIDE PXD038203), staged (scripts/00_stage_study_sites.py). **Terms:** CC0 (PRIDE PXD038203 Table S1), the article is not open access.

> MaxQuant ADMA-antibody sites with a 31-residue window, localisation probability >= 0.75, 1 % FDR. Tumours grown in mice: a peptide identical in mouse and human can come from mouse stroma (flag xenograft_only).

403 records: 403 placed on 403 sites (359 of them in the previous version, 46 reported by no other source). Every record is placed.

```text
record 66258: reported A0AV96 at R394 (localisation probability 1, 1 source row, tables Table S1 Identified ADMA sites, patient-derived xenograft in mouse)
placed on A0AV96:R394 (sequence_unique_canonical; evidence check: agrees); the window of 31 residues occurs once in the reviewed proteome
evidence              RGRGRGAAGNRAPGP[R]GSYLGGYSAGRGIYS
reference   373 VKAGSIRGRGRGAAGNRAPGP[R]GSYLGGYSAGRGIYSRYHEGK 415
```
