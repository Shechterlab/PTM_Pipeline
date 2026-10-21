"""Disease families for the DisGeNET enrichment (used by scripts/14_disease_enrichment.py and scripts/34).

Each DisGeNET term is assigned to at most one of eight audited families by the first rule its name matches; a few names
are assigned by hand (CURATED_OVERRIDES). The rules are those of the original pipeline (scripts/34; the disease analysis is
not part of the paper), except that v1.1.0 adds TDP-43 proteinopathy, Pick disease and anterior horn cell disease to the
ALS / motor neuron / FTD rule. A term no audited rule matches goes to the first of the additional families
(ADDITIONAL_FAMILY_RULES, rule-based and not audited) whose rule it matches; those are tested apart from the audited
families. Terms that match no rule are "Other / unassigned" and are not tested.
"""
from __future__ import annotations

import re

# First match wins. This precedence is the reproducible primary assignment rule.
FAMILY_RULES = [
    # v1.1.0 adds TDP-43 proteinopathy, Pick disease (not Niemann-Pick) and anterior horn cell disease, which the
    # published rule missed.
    ("ALS / motor neuron / FTD", r"amyotrophic|\bals\b|motor neuron|frontotemporal|\bftd\b|bulbar palsy"
                                 r"|tdp-?43 proteinopath|(?<!niemann-)pick(?:'s)? disease|anterior horn"),
    ("Alzheimer / tau / dementia", r"alzheimer|tauopath|\bdementia\b|cognitive disorder"),
    ("Parkinsonism", r"parkinson"),
    ("Ataxia / cerebellar degeneration", r"ataxia|spinocerebellar|friedreich|cerebellar"),
    ("Myopathy / muscle / cardiomyopathy", r"myopath|muscular dystrophy|cardiomyopath|rhabdomyolysis|myositis"),
    ("Cancer / neoplasia", r"cancer|carcinoma|sarcoma|glioma|leukemia|lymphoma|melanoma|tumou?r|adenoma|neoplas|myeloma"),
    ("Immune / inflammatory", r"autoimmune|inflamm|lupus|arthritis|colitis|crohn|asthma|dermatitis|psoriasis"),
    ("Developmental / congenital syndromes", r"congenital|dysmorph|cranio|palate|facies|syndrom|hypoplasia|epicanthus|forehead|midface|feeding difficulties"),
]

# Applied only to terms no audited rule (or override) assigns; first match wins. Mostly clinical signs (HPO terms).
ADDITIONAL_FAMILY_RULES = [
    ("Cancer-related (other)", r"leukaemia|malignan|teratoma|astrocytoma|blastoma|neurofibromatos|fibromatosis|lymphangioleiomyomatosis"
                               r"|microsatellite instability|replication error|somatic mutation|myelodysplas|myeloproliferative|blast phase"
                               r"|excess blasts"),
    ("Neurodegeneration and neuromuscular (other)", r"neurodegenerat|proteinopath|pick disease|huntington|chorea|leukodystroph"
                               r"|leukoencephalopath|spastic paraplegia|muscular atrophy|anterior horn|machado-joseph|gliosis|brain atrophy"
                               r"|cortical atrophy|alexander disease|leigh|motor and sensory|charcot-marie|neuropath|amyotroph|fasciculation"
                               r"|bulbar|nerve conduction|denervation|neuromuscular|myoton|muscle weakness|muscular weakness|weakness of"
                               r"|muscle fiber|muscle atrophy|dystroph|rimmed vacuoles|myoglobinuria|creatine phosphokinase|onion bulb|little's disease"),
    ("Neurodevelopmental and behavioural", r"intellectual|mental retardation|mental deficiency|autis|speech|language|learning|developmental delay"
                               r"|development delay|behavio|school|aggressive|mood|agitation|obsessive|self-injurious|stereotyp|mutism"
                               r"|eye contact|psychomotor|microcephal|macrocephal|large head|polymicrogyria|lissenceph|pachygyria"
                               r"|corpus callosum|gyral|myelination|electroencephalogram|ability to walk|head control|dysphasia"),
    ("Neurological signs", r"reflex|areflexia|babinski|pyramidal|spastic|spasticity|paresis|paraly|hypertonia|hypotonia|spasm|tremor"
                           r"|gait|walking|falls|ataxic|coordination|ophthalmoplegia|pursuit|drooling|neurologic|cranial nerve"
                           r"|syringomyelia|arnold chiari|arachnoid|sensory impairment|pain sensation|dysphonia|voice|cry\b"),
    ("Cardiovascular", r"cardi|heart|valv|aort|arrhythm|tachycardia|atrioventricular|septal|fallot|foramen ovale|pulmonary stenosis"
                       r"|coarctation|aneurysm|arteriovenous|capillary malformation|port-wine|lymphedema|endocardi"),
    ("Blood and haematological", r"anemi|anaemi|platelet|bruis|corpuscular|reticulocyte|marrow|spleen|pallor|purpura|bleeding"),
    ("Metabolic and mitochondrial", r"metabol|mitochondri|acidosis|lactate|pyruvate|ammon|phenylketon|glycogen|enzymopathy|amino acid"
                                    r"|protein deficiency|insulin|obesity|weight|calcitonin"),
    ("Eye and ear", r"hearing|cataract|ophthalm|\beye|eyes\b|visual|esotropia|exotropia|sclera|hyperopia|ptosis|blepharo|chorioretin"
                    r"|extraocular|exophthalmos|phthisis|embryotoxon|auditory|orbit"),
    ("Kidney, urinary and reproductive", r"kidney|renal|urinary|ureter|urethr|genital|labia|menarche|menopause|puberty|metrorrhagia|nipple"),
    ("Digestive and liver", r"liver|hepat|ulcer|stomach|pyloric|intestin|volvulus|\banus\b|\banal\b|reflux|gingiv|oral cavity|hernia"),
    ("Infection", r"infection|virus|viral|\bhiv\b|\baids\b|parvovirus|mycoplasma|adenovirus|enterovirus|coinfection"),
    ("Skin, hair and nails", r"skin|epiderm|keratosis|xeroderma|cafe-au-lait|hirsut|hypertrichosis|cutis|dermatoglyph|pachyonychia|nail"
                             r"|\bhair|eyelash|eyebrow|synophrys"),
    ("Respiratory", r"respirat|dyspnea|cough|lung|asbestosis|laryng|trache|vocal cord|asphyxia"),
    ("Pregnancy, birth and growth", r"prematur|gestational|birth|fetal|polyhydramnios|nuchal|growth retardation|failure to thrive"
                                    r"|intrauterine|overgrowth|progeria|suck|short stature|proportionate|slender build|percentile"),
    ("Craniofacial, skeletal and limb features", r"face|facial|philtrum|\blip|lips\b|nasal|nose|\bear|ears\b|earlobe|pinna|helix|macrotia"
                                    r"|palpebral|forehead|chin|mouth|microstomia|macrostomia|dental|dentition|tooth|teeth|hypodontia|malocclusion"
                                    r"|bite\b|retrognathia|mandib|hairline|cephal|occiput|cheek|telecanthus|spine|spina|kyphos|lordosis|scolio"
                                    r"|joint|contracture|finger|\btoe|toes\b|toenail|hand\b|hands\b|foot\b|feet\b|flatfoot|clinodactyly"
                                    r"|syndactyly|arachnodactyly|talipes|thumb|phalan|palmar|micromelia|hip dysplasia|scapula|shoulder|thorax"
                                    r"|neck|head\b|malar|sacral|sprengel|ankle|elbow|bone|osteo|epiphys|metaphys|synostosis|arthrogryposis|sandal gap|meniscal|vermilion|eversion"),
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
        extra = [fam for fam, pat in ADDITIONAL_FAMILY_RULES if re.search(pat, name, re.I)]
        if extra:
            return extra[0], "included", "additional_rule", f"No audited family matched; first matching additional family ({extra[0]}).", " | ".join(extra[1:])
        return "Other / unassigned", "unassigned", "no_match", "No overview-family label rule matched.", ""
    return matches[0], "included", "precedence_rule", f"Primary assignment is the first matching family in precedence order ({matches[0]}).", " | ".join(matches[1:])


AUDITED_FAMILIES = [fam for fam, _ in FAMILY_RULES]
ADDITIONAL_FAMILIES = [fam for fam, _ in ADDITIONAL_FAMILY_RULES]
ALL_FAMILIES = AUDITED_FAMILIES + ADDITIONAL_FAMILIES
