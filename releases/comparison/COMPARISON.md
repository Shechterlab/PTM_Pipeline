# Analyses of v1.0.0, v1.1.0, v1.2.0, v1.2.0-lc compared

Each release's sites were analysed with the same code (`hpc/release_analyses.sh`; the commit is in each `analyses/README.md`), so the differences below come from the sites, not the methods. A test is called enriched or depleted at BH q <= 0.05 (the effect is the log2 fold enrichment, or the log2 odds ratio). `comparison_v1.0.0_v1.1.0_v1.2.0_v1.2.0-lc.tsv` has every test side by side.

## Headline numbers

| | v1.0.0 | v1.1.0 | v1.2.0 | v1.2.0-lc |
|---|---:|---:|---:|---:|
| sites analysed | 12,495 | 13,365 | 11,601 | 13,445 |
| in a disordered region | 48.5% | 49.5% | 50.6% | 49.6% |
| at a disorder boundary | 23.6% | 23.1% | 22.6% | 23.1% |
| inside a domain | 34.1% | 33.6% | 32.1% | 33.4% |
| RG at the site | 27.7% | 27.8% | 30.2% | 27.9% |
| RG, odds ratio vs all arginines | 5.89 | 5.97 | 6.68 | 5.99 |
| RGG at the site | 7.0% | 7.0% | 7.8% | 7.0% |
| RGG, odds ratio vs all arginines | 18.90 | 19.56 | 21.42 | 19.28 |
| GAR at the site | 0.7% | 0.7% | 0.8% | 0.8% |
| GAR, odds ratio vs all arginines | 1.56 | 1.65 | 1.87 | 1.74 |
| RXR at the site | 16.4% | 17.5% | 16.7% | 17.5% |
| RXR, odds ratio vs all arginines | 1.13 | 1.22 | 1.15 | 1.23 |
| CARM1_proline_rich at the site | 50.1% | 49.6% | 50.2% | 50.0% |
| CARM1_proline_rich, odds ratio vs all arginines | 1.40 | 1.38 | 1.41 | 1.40 |
| another site within 5 aa | 27.2% | 30.3% | 23.8% | 30.3% |
| another site within 5 aa, obs/null | 2.16 | 2.19 | 1.84 | 2.19 |
| another site within 10 aa | 41.7% | 44.8% | 37.8% | 44.7% |
| another site within 10 aa, obs/null | 1.93 | 1.94 | 1.78 | 1.94 |
| another site within 20 aa | 51.4% | 54.3% | 48.2% | 54.2% |
| another site within 20 aa, obs/null | 1.63 | 1.63 | 1.57 | 1.64 |
| proteins with a dense window (any width) | 180 | 259 | 217 | 270 |
| proteins with a dense 10-aa window | 34 | 58 | 45 | 57 |
| proteins with a dense 20-aa window | 58 | 100 | 75 | 100 |
| proteins with a dense 50-aa window | 70 | 116 | 102 | 125 |
| methylated proteins in a condensate | 45.2% | 44.7% | 45.5% | 44.6% |

## Tests by analysis

| Analysis | Tested | Enriched v1.0.0 | Enriched v1.1.0 | Enriched v1.2.0 | Enriched v1.2.0-lc | Same call in all | Spearman r of effects (first vs last) |
|---|---:|---:|---:|---:|---:|---:|---:|
| functional class | 13 | 10 | 10 | 10 | 10 | 13 | 0.99 |
| PANTHER group | 11 | 4 | 4 | 5 | 5 | 10 | 1.00 |
| PANTHER class | 25 | 7 | 8 | 8 | 9 | 22 | 1.00 |
| motif family at the site | 11 | 7 | 8 | 7 | 8 | 10 | 1.00 |
| disorder context | 4 | 1 | 1 | 1 | 1 | 3 | 0.40 |
| domain context | 5 | 3 | 3 | 3 | 3 | 5 | 1.00 |
| nearest domain class | 14 | 6 | 7 | 7 | 7 | 13 | 0.98 |
| CD-CODE condensate | 40 | 10 | 10 | 13 | 10 | 37 | 0.95 |
| HPA location | 48 | 11 | 11 | 9 | 7 | 43 | 0.97 |
| disease family (audited) | 8 | 7 | 7 | 8 | 7 | 7 | 0.98 |
| functional subclass | 27 | 22 | 22 | 22 | 22 | 27 | 1.00 |
| PANTHER subclass | 129 | 32 | 33 | 31 | 35 | 120 | 0.99 |
| disease family (additional) | 15 | 15 | 15 | 15 | 15 | 15 | 1.00 |

### Functional class

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| Translation / ribosome / biogenesis | **1.15** (3.8e-72) | **1.15** (8.6e-76) | **1.15** (6.8e-69) | **1.13** (3.9e-72) |
| Nuclear scaffold / pore / lamina | **0.989** (1.6e-14) | **0.996** (2.1e-15) | **0.928** (4.9e-12) | **1** (1.1e-15) |
| RNA-binding / RNP | **0.893** (8.9e-124) | **0.883** (1.7e-125) | **0.901** (6.2e-121) | **0.884** (7.2e-127) |
| Cell cycle / centrosome / spindle | **0.751** (1.7e-07) | **0.774** (3.1e-08) | **0.816** (1.5e-08) | **0.766** (4e-08) |
| DNA repair / genome maintenance | **0.734** (1.8e-11) | **0.72** (2.1e-11) | **0.763** (7e-12) | **0.741** (3.9e-12) |
| Proteostasis / ubiquitin / chaperone | **0.651** (1.1e-23) | **0.651** (1.4e-24) | **0.666** (6.5e-24) | **0.653** (5.8e-25) |
| Cytoplasmic / trafficking / cytoskeleton | **0.517** (2.6e-29) | **0.523** (2.4e-31) | **0.512** (1.3e-27) | **0.51** (8.7e-30) |
| Mitochondrial / metabolism | **0.445** (7e-12) | **0.438** (5.9e-12) | **0.447** (1.4e-11) | **0.434** (7.4e-12) |
| Chromatin / transcription | **0.349** (1.2e-20) | **0.35** (1.2e-21) | **0.353** (2.6e-20) | **0.362** (2.7e-23) |
| Signaling / enzymes | **0.365** (6e-37) | **0.367** (7.1e-39) | **0.373** (3e-37) | **0.361** (9.1e-38) |
| Other / unclassified | -0.366 (1) | -0.369 (1) | -0.34 (1) | -0.36 (1) |
| Membrane / receptor / extracellular | -0.891 (1) | -0.886 (1) | -0.889 (1) | -0.887 (1) |
| Poorly characterized / specialized | -0.94 (1) | -0.916 (1) | -1.05 (1) | -0.958 (1) |

### PANTHER group

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| RNA processing and translation | **1.28** (2.6e-129) | **1.27** (2.1e-132) | **1.29** (4.3e-128) | **1.26** (4e-131) |
| Scaffolds, adaptors and chaperones | **0.472** (6.2e-10) | **0.456** (1.2e-09) | **0.425** (9.5e-08) | **0.444** (3.6e-09) |
| Transcription, chromatin and DNA | **0.196** (0.00032) | **0.204** (0.00012) | **0.169** (0.0031) | **0.221** (2.3e-05) |
| Cytoskeleton, adhesion and matrix | **0.146** (0.046) | **0.145** (0.044) | **0.159** (0.025) | **0.166** (0.018) |
| Protein-modifying enzymes † | 0.113 (0.078) | 0.12 (0.056) | **0.183** (0.0047) | **0.13** (0.037) |
| Metabolic enzymes | 0.033 (0.53) | 0.0337 (0.52) | 0.0252 (0.63) | 0.0444 (0.4) |
| No PANTHER protein class | -0.137 (1) | -0.142 (1) | -0.125 (1) | -0.137 (1) |
| Transport and membrane traffic | -0.344 (1) | -0.332 (1) | -0.396 (1) | -0.385 (1) |
| Other PANTHER classes | -0.453 (1) | -0.497 (1) | -0.405 (1) | -0.505 (1) |
| Signalling: receptors, ligands and modulators | -0.943 (1) | -0.918 (1) | -0.978 (1) | -0.941 (1) |
| Immunity and defence | -2.6 (1) | -2.58 (1) | -2.75 (1) | -2.65 (1) |

### PANTHER class

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| translational protein | **1.41** (1.9e-48) | **1.42** (2.2e-51) | **1.38** (5.4e-43) | **1.38** (4.9e-48) |
| RNA metabolism protein | **1.22** (2.4e-78) | **1.2** (8.5e-79) | **1.26** (6.2e-82) | **1.2** (2e-80) |
| chromatin/chromatin-binding, or -regulatory protein | **1.08** (3e-22) | **1.09** (6.5e-24) | **0.946** (2.4e-15) | **1.09** (1.2e-23) |
| chaperone | **1.09** (1e-15) | **1.06** (2.8e-15) | **1.05** (1e-13) | **1.03** (4e-14) |
| cytoskeletal protein | **0.75** (4.3e-18) | **0.743** (2.3e-18) | **0.749** (3.8e-17) | **0.752** (4.7e-19) |
| DNA metabolism protein | **0.664** (2e-05) | **0.638** (3.4e-05) | **0.657** (4e-05) | **0.63** (4.2e-05) |
| membrane traffic protein † | 0.25 (0.064) | **0.291** (0.021) | **0.31** (0.018) | **0.307** (0.015) |
| scaffold/adaptor protein † | **0.232** (0.02) | **0.225** (0.021) | 0.182 (0.077) | **0.224** (0.019) |
| protein modifying enzyme † | 0.113 (0.099) | 0.12 (0.071) | **0.183** (0.0062) | **0.13** (0.047) |
| protein-binding activity modulator | 0.0853 (0.46) | 0.108 (0.29) | 0.0553 (0.72) | 0.0931 (0.39) |
| metabolite interconversion enzyme | 0.033 (0.66) | 0.0337 (0.64) | 0.0252 (0.78) | 0.0444 (0.5) |
| calcium-binding protein | -0.0334 (1) | -0.0205 (1) | 0.0146 (1) | -0.0284 (1) |
| extracellular matrix protein | -0.273 (1) | -0.262 (1) | -0.225 (1) | -0.118 (1) |
| gene-specific transcriptional regulator | -0.166 (1) | -0.153 (1) | -0.15 (1) | -0.118 (1) |
| No PANTHER protein class | -0.137 (1) | -0.142 (1) | -0.125 (1) | -0.137 (1) |
| viral or transposable element protein | -0.402 (1) | -0.446 (1) | -0.354 (1) | -0.454 (1) |
| cell junction protein | -0.581 (1) | -0.473 (1) | -0.533 (1) | -0.481 (1) |
| transfer/carrier protein | -0.539 (1) | -0.583 (1) | -0.636 (1) | -0.662 (1) |
| transporter | -0.657 (1) | -0.662 (1) | -0.798 (1) | -0.769 (1) |
| cell adhesion molecule | -0.836 (1) | -0.811 (1) | -0.788 (1) | -0.785 (1) |
| storage protein | -0.844 (1) | -0.888 (1) | -0.796 (1) | -0.896 (1) |
| structural protein | -1.09 (1) | -1.14 (1) | -1.04 (1) | -1.14 (1) |
| intercellular signal molecule | -1.99 (1) | -1.97 (1) | -1.94 (1) | -2.04 (1) |
| transmembrane signal receptor | -2.39 (1) | -2.34 (1) | -2.53 (1) | -2.38 (1) |
| defense/immunity protein | -2.6 (1) | -2.58 (1) | -2.75 (1) | -2.65 (1) |

### Motif family at the site

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| RGG | **4.23** (0) | **4.36** (0) | **4.47** (0) | **4.33** (0) |
| GRG | **3.88** (0) | **4.1** (0) | **4.22** (0) | **4.1** (0) |
| RG | **2.6** (0) | **2.63** (0) | **2.79** (0) | **2.64** (0) |
| DR | **1.52** (1.6e-210) | **1.47** (3.3e-207) | **1.56** (2e-212) | **1.48** (3.5e-211) |
| GAR | **0.63** (0.00028) | **0.707** (2e-05) | **0.876** (2.4e-07) | **0.791** (9.9e-07) |
| CARM1_proline_rich | **0.504** (2.3e-79) | **0.477** (1.1e-75) | **0.507** (9.3e-75) | **0.502** (2.4e-84) |
| CARM1_hydrophobic_proline | **0.472** (6.3e-68) | **0.431** (1.4e-60) | **0.479** (2.4e-65) | **0.459** (7e-69) |
| RXR † | 0.0337 (0.35) | **0.147** (2e-05) | 0.0637 (0.089) | **0.152** (9.6e-06) |
| ER | **-0.717** (5e-37) | **-0.728** (2.3e-40) | **-0.778** (1.5e-39) | **-0.779** (1e-45) |
| RD | **-0.842** (6.4e-33) | **-0.863** (2.6e-36) | **-0.88** (9.7e-33) | **-0.863** (1.9e-36) |
| RE | **-1.31** (2.1e-94) | **-1.25** (2.5e-93) | **-1.47** (9.5e-104) | **-1.32** (2.5e-102) |

### Disorder context

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| disordered | **0.934** (1.3e-261) | **1.01** (0) | **1.04** (4.4e-306) | **1.01** (0) |
| disorder_boundary | **-0.0788** (0.016) | **-0.125** (5.3e-05) | **-0.154** (3.2e-06) | **-0.117** (0.00016) |
| ordered | **-0.928** (1e-232) | **-0.979** (7.6e-275) | **-1** (6.7e-250) | **-0.979** (2e-275) |
| no_disorder_annotation † | 0.601 (0.45) | 0.546 (0.46) | 0.922 (0.13) | **-1.01** (0.022) |

### Domain context

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| no_domain_annotation | **0.515** (9.7e-37) | **0.554** (4.1e-45) | **0.564** (2.7e-40) | **0.56** (4.6e-46) |
| distal | **0.498** (1.2e-64) | **0.538** (2.2e-80) | **0.605** (2.8e-90) | **0.542** (2e-82) |
| boundary | **0.151** (0.0014) | **0.148** (0.0012) | **0.16** (0.0011) | **0.135** (0.0032) |
| inter_domain_linker | **-0.332** (6.4e-16) | **-0.345** (4.4e-18) | **-0.33** (7.8e-15) | **-0.351** (7.5e-19) |
| in_domain | **-0.552** (2.9e-88) | **-0.598** (3.8e-110) | **-0.681** (8.5e-123) | **-0.597** (3.7e-110) |

### Nearest domain class

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| Sm / Lsm | **1.99** (0.00065) | **2.37** (1e-05) | **2.09** (0.00039) | **2.39** (9e-06) |
| RRM / KH RNA-binding | **1.81** (1.8e-210) | **1.85** (4.6e-230) | **1.86** (1.7e-209) | **1.87** (3.9e-235) |
| Ribosome / RNP structural | **1.64** (6.9e-39) | **1.78** (1.6e-50) | **1.75** (1.2e-43) | **1.82** (1.7e-52) |
| dsRNA-binding | **0.8** (0.042) | **1.05** (0.0038) | **1.2** (0.0012) | **1.06** (0.0035) |
| No domain annotation | **0.515** (2e-36) | **0.554** (8.6e-45) | **0.564** (5.6e-40) | **0.56** (9.7e-46) |
| RNA-associated zinc finger † | 0.202 (0.23) | **0.396** (0.0066) | **0.376** (0.013) | **0.373** (0.011) |
| DNA-binding | **0.301** (0.024) | **0.352** (0.0039) | **0.378** (0.003) | **0.303** (0.013) |
| Chromatin catalytic / remodeler | 0.0548 (0.75) | 0.0621 (0.68) | -0.0602 (0.74) | 0.0633 (0.68) |
| RNA helicase | -0.115 (0.25) | -0.152 (0.11) | -0.179 (0.089) | -0.152 (0.12) |
| Catalytic enzyme domain | **-0.301** (1.3e-08) | **-0.36** (1.4e-12) | **-0.352** (1.1e-10) | **-0.323** (2.1e-10) |
| Protein-protein repeat / scaffold | **-0.387** (2.1e-07) | **-0.385** (1.2e-07) | **-0.393** (4.2e-07) | **-0.362** (4.3e-07) |
| Other domain | **-0.421** (1.8e-53) | **-0.425** (3.2e-58) | **-0.441** (6.8e-55) | **-0.432** (1.6e-60) |
| Repeat / superfamily | **-0.36** (1.1e-07) | **-0.435** (9.3e-12) | **-0.397** (7.6e-09) | **-0.457** (5.4e-13) |
| Chromatin reader | **-0.563** (3.4e-05) | **-0.55** (2.6e-05) | **-0.522** (0.00018) | **-0.569** (1.1e-05) |

### CD-CODE condensate

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| 3190BB68 † | 0.957 (0.057) | 0.93 (0.064) | **0.997** (0.048) | 0.926 (0.065) |
| 373365D7 † | 0.957 (0.057) | 0.93 (0.064) | **0.997** (0.048) | 0.926 (0.065) |
| BEB5E62D † | 0.957 (0.057) | 0.93 (0.064) | **0.997** (0.048) | 0.926 (0.065) |
| A40D6D0A | **0.832** (0.021) | **0.93** (0.0025) | **0.997** (0.0013) | **0.926** (0.0025) |
| F85A2E29 | **0.839** (6.2e-07) | **0.811** (1.1e-06) | **0.919** (2.5e-08) | **0.848** (1.1e-07) |
| D8E9712B | **0.765** (0.017) | **0.837** (0.0025) | **0.904** (0.0013) | **0.833** (0.0026) |
| 62EA6512 | **0.751** (0.022) | **0.723** (0.028) | **0.791** (0.015) | **0.72** (0.026) |
| 232F8A39 | **0.467** (6.4e-34) | **0.455** (2.3e-33) | **0.483** (1.4e-34) | **0.459** (2.6e-34) |
| DEE660B4 | **0.41** (8.4e-41) | **0.395** (3.8e-39) | **0.428** (2.2e-42) | **0.397** (7.7e-40) |
| 1A18FFC4 | **0.383** (0.022) | **0.355** (0.035) | **0.453** (0.0054) | **0.382** (0.02) |
| 804901D1 | **0.351** (0.00011) | **0.345** (0.0001) | **0.433** (8.1e-07) | **0.362** (3.3e-05) |
| B5B9A610 | **0.354** (2e-09) | **0.346** (2.1e-09) | **0.336** (4.2e-08) | **0.343** (2.9e-09) |
| 91857CE7 | **0.297** (6.1e-24) | **0.303** (5.5e-26) | **0.293** (8e-22) | **0.3** (1.8e-25) |

(Only the categories called in some release are listed; 27 others are not significant in any.)

### HPA location

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| Endosomes † | **1.07** (0.041) | **1.03** (0.041) | 0.958 (0.08) | 1.02 (0.053) |
| Mitotic chromosome | **0.745** (0.0049) | **0.785** (0.0014) | **0.829** (0.00098) | **0.776** (0.0016) |
| Primary cilium tip | **0.702** (0.0093) | **0.662** (0.014) | **0.607** (0.033) | **0.609** (0.035) |
| Nucleoli rim | **0.515** (0.0093) | **0.572** (0.002) | **0.558** (0.0051) | **0.54** (0.0046) |
| Calyx † | **0.585** (0.041) | **0.593** (0.039) | **0.629** (0.029) | 0.537 (0.06) |
| Centriolar satellite † | **0.385** (0.041) | **0.369** (0.048) | 0.379 (0.053) | 0.336 (0.077) |
| Basal body † | **0.308** (0.041) | **0.295** (0.039) | 0.281 (0.053) | 0.259 (0.072) |
| Nuclear speckles | **0.233** (0.041) | **0.23** (0.041) | **0.267** (0.029) | **0.248** (0.037) |
| Cytosol | **0.215** (2.5e-16) | **0.204** (2.9e-15) | **0.208** (1.1e-14) | **0.202** (5.2e-15) |
| Nuclear bodies † | 0.165 (0.11) | 0.164 (0.11) | **0.224** (0.033) | 0.185 (0.072) |
| Nucleoli | **0.164** (0.041) | **0.166** (0.039) | **0.189** (0.028) | **0.166** (0.04) |
| Nucleoplasm | **0.134** (7.7e-09) | **0.134** (3.8e-09) | **0.138** (7.1e-09) | **0.138** (9.6e-10) |

(Only the categories called in some release are listed; 36 others are not significant in any.)

### Disease family (audited)

| Category | v1.0.0: log2 effect (q) | v1.1.0: log2 effect (q) | v1.2.0: log2 effect (q) | v1.2.0-lc: log2 effect (q) |
|---|---|---|---|---|
| Ataxia / cerebellar degeneration | **0.447** (6.1e-20) | **0.45** (5.3e-21) | **0.48** (2e-22) | **0.44** (3.2e-20) |
| ALS / motor neuron / FTD | **0.373** (1.5e-11) | **0.365** (1.4e-11) | **0.381** (9.6e-12) | **0.341** (3e-10) |
| Myopathy / muscle / cardiomyopathy | **0.268** (1.5e-11) | **0.264** (1.4e-11) | **0.283** (2.5e-12) | **0.258** (4.3e-11) |
| Alzheimer / tau / dementia | **0.196** (3.4e-08) | **0.191** (4.3e-08) | **0.201** (2.6e-08) | **0.175** (4.6e-07) |
| Parkinsonism | **0.204** (1.3e-05) | **0.188** (4.5e-05) | **0.186** (9.2e-05) | **0.161** (0.00041) |
| Developmental / congenital syndromes | **0.148** (9.4e-17) | **0.141** (9.3e-16) | **0.147** (5.9e-16) | **0.135** (1.1e-14) |
| Cancer / neoplasia | **0.102** (2.6e-29) | **0.0992** (8.3e-29) | **0.101** (4.6e-28) | **0.0981** (2e-28) |
| Immune / inflammatory † | 0.0387 (0.054) | 0.0374 (0.057) | **0.0466** (0.028) | 0.0309 (0.096) |

† the call (enriched, depleted, not significant) differs between releases. Bold: q <= 0.05.
