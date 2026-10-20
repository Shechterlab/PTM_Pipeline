"""Disease families for the DisGeNET enrichment (used by scripts/14_disease_enrichment.py and scripts/34).

Each DisGeNET term is assigned to at most one of eight families by the first rule its name matches; a few names are
assigned by hand (CURATED_OVERRIDES). Terms that match no rule are "Other / unassigned" and are not tested.
"""
from __future__ import annotations

import re

# First match wins. This precedence is the reproducible primary assignment rule.
FAMILY_RULES = [
    ("ALS / motor neuron / FTD", r"amyotrophic|\bals\b|motor neuron|frontotemporal|\bftd\b|bulbar palsy"),
    ("Alzheimer / tau / dementia", r"alzheimer|tauopath|\bdementia\b|cognitive disorder"),
    ("Parkinsonism", r"parkinson"),
    ("Ataxia / cerebellar degeneration", r"ataxia|spinocerebellar|friedreich|cerebellar"),
    ("Myopathy / muscle / cardiomyopathy", r"myopath|muscular dystrophy|cardiomyopath|rhabdomyolysis|myositis"),
    ("Cancer / neoplasia", r"cancer|carcinoma|sarcoma|glioma|leukemia|lymphoma|melanoma|tumou?r|adenoma|neoplas|myeloma"),
    ("Immune / inflammatory", r"autoimmune|inflamm|lupus|arthritis|colitis|crohn|asthma|dermatitis|psoriasis"),
    ("Developmental / congenital syndromes", r"congenital|dysmorph|cranio|palate|facies|syndrom|hypoplasia|epicanthus|forehead|midface|feeding difficulties"),
]

# Exact-name overrides are deliberately separate from label matching.
CURATED_OVERRIDES = {
    "Ataxia Telangiectasia": (
        "Excluded: genome-instability syndrome",
        "excluded_from_overview",
        "Curated override: excluded rather than counted as general ataxia; genome-instability family is outside the 8-family overview.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE AND FRONTOTEMPORAL DEMENTIA": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "Inclusion Body Myopathy with Early-Onset Paget Disease with or without Frontotemporal Dementia 1": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE WITH OR WITHOUT FRONTOTEMPORAL DEMENTIA 2": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE WITH OR WITHOUT FRONTOTEMPORAL DEMENTIA 3": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
}

FOCUS_TERMS = [
    ("C4024896", "Motor neuron atrophy"),
    ("C0085084", "Motor neuron disease"),
    ("C0338451", "Frontotemporal dementia (FTD)"),
    ("C4551993", "Familial ALS"),
    ("C0949664", "Tauopathies"),
    ("C0751072", "Frontotemporal lobar degeneration"),
]



def classify(name: str) -> tuple[str, str, str, str, str]:
    """(family, audit status, method, reason, other matching families) of a DisGeNET term name."""
    if name in CURATED_OVERRIDES:
        fam, status, reason = CURATED_OVERRIDES[name]
        return fam, status, "curated_override", reason, ""
    matches = [fam for fam, pat in FAMILY_RULES if re.search(pat, name, re.I)]
    if not matches:
        return "Other / unassigned", "unassigned", "no_match", "No overview-family label rule matched.", ""
    return matches[0], "included", "precedence_rule", f"Primary assignment is the first matching family in precedence order ({matches[0]}).", " | ".join(matches[1:])
