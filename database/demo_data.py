"""
Curated demo dataset for OpenTargetAI.

A small, scientifically plausible set of well-characterised drugs, their primary
protein targets and approximate potencies (drawn from public knowledge of
ChEMBL/BindingDB bioactivities). It is intentionally organised into pharmacological
families with several analogues per target so that — without downloading the full
multi-gigabyte ChEMBL/BindingDB exports — similarity search, target prediction,
the explainer and the benchmark all return meaningful results out of the box.

This is *demonstration* data: potency values are rounded, representative figures
and must not be used for real decision making.
"""

from __future__ import annotations

# (name, SMILES, ChEMBL id or None)
COMPOUNDS: list[tuple[str, str, str | None]] = [
    # NSAIDs / coxibs
    ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "CHEMBL521"),
    ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1", "CHEMBL154"),
    ("Ketoprofen", "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1", "CHEMBL558"),
    ("Flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1", "CHEMBL563"),
    ("Diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "CHEMBL139"),
    ("Indomethacin", "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1", "CHEMBL6"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", "CHEMBL25"),
    ("Celecoxib", "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1", "CHEMBL118"),
    ("Mefenamic acid", "Cc1cccc(C)c1Nc1ccccc1C(=O)O", "CHEMBL536"),
    # Statins
    (
        "Atorvastatin",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        "CHEMBL1487",
    ),
    ("Simvastatin", "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12", "CHEMBL1064"),
    (
        "Pravastatin",
        "CCC(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC(O)CC(O)CC(=O)O)C21",
        "CHEMBL1144",
    ),
    (
        "Rosuvastatin",
        "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1C=CC(O)CC(O)CC(=O)O",
        "CHEMBL1496",
    ),
    (
        "Fluvastatin",
        "CC(C)n1c(C=CC(O)CC(O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21",
        "CHEMBL1164",
    ),
    # SSRIs
    ("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "CHEMBL41"),
    ("Paroxetine", "Fc1ccc(C2CCNCC2COc2ccc3c(c2)OCO3)cc1", "CHEMBL490"),
    ("Sertraline", "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21", "CHEMBL809"),
    ("Citalopram", "CN(C)CCCC1(c2ccc(C#N)cc2)OCc2cc(F)ccc21", "CHEMBL549"),
    ("Fluvoxamine", "COCCCCC(=NOCCN)c1ccc(C(F)(F)F)cc1", "CHEMBL1407"),
    # Beta blockers
    ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", "CHEMBL27"),
    ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", "CHEMBL24"),
    ("Metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1", "CHEMBL13"),
    ("Bisoprolol", "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1", "CHEMBL526"),
    # ACE inhibitors
    ("Lisinopril", "NCCCCC(NC(CCc1ccccc1)C(=O)O)C(=O)N1CCCC1C(=O)O", "CHEMBL1237"),
    ("Captopril", "CC(CS)C(=O)N1CCCC1C(=O)O", "CHEMBL1208"),
    ("Enalapril", "CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(=O)O", "CHEMBL578"),
    # Kinase inhibitors
    (
        "Imatinib",
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "CHEMBL941",
    ),
    (
        "Dasatinib",
        "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1",
        "CHEMBL1421",
    ),
    ("Gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", "CHEMBL939"),
    ("Erlotinib", "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", "CHEMBL553"),
    # Estrogen receptor modulators
    ("Tamoxifen", "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1", "CHEMBL83"),
    (
        "Raloxifene",
        "O=C(c1ccc(O)cc1)c1c(-c2ccc(OCCN3CCCCC3)cc2)sc2cc(O)ccc12",
        "CHEMBL81",
    ),
    # PDE5 inhibitors
    (
        "Sildenafil",
        "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
        "CHEMBL192",
    ),
    ("Tadalafil", "O=C1N(C)CC(=O)N2C1Cc1c2cc2c(c1)OCO2", "CHEMBL779"),
    # H1 antihistamines
    ("Diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", "CHEMBL657"),
    ("Loratadine", "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1", "CHEMBL998"),
    ("Cetirizine", "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1", "CHEMBL1000"),
]

# (UniProt, gene, name, organism, target_class, description, disease_associations)
TARGETS: list[tuple[str, str, str, str, str, str, str]] = [
    (
        "P23219",
        "PTGS1",
        "Prostaglandin G/H synthase 1 (COX-1)",
        "Homo sapiens",
        "Enzyme (Oxidoreductase)",
        "Constitutive cyclooxygenase catalysing the committed step in "
        "prostaglandin biosynthesis; the principal target of classical NSAIDs.",
        "Pain; Inflammation; Cardiovascular (gastric side-effects)",
    ),
    (
        "P35354",
        "PTGS2",
        "Prostaglandin G/H synthase 2 (COX-2)",
        "Homo sapiens",
        "Enzyme (Oxidoreductase)",
        "Inducible cyclooxygenase up-regulated during inflammation; the target "
        "of selective coxib anti-inflammatory drugs.",
        "Inflammation; Pain; Colorectal cancer",
    ),
    (
        "P04035",
        "HMGCR",
        "3-hydroxy-3-methylglutaryl-CoA reductase",
        "Homo sapiens",
        "Enzyme (Oxidoreductase)",
        "Rate-limiting enzyme of cholesterol biosynthesis; the target of statins.",
        "Hypercholesterolaemia; Atherosclerosis; Cardiovascular disease",
    ),
    (
        "P31645",
        "SLC6A4",
        "Sodium-dependent serotonin transporter (SERT)",
        "Homo sapiens",
        "Transporter",
        "Reuptake transporter terminating serotonergic signalling; the target of "
        "SSRIs and tricyclic antidepressants.",
        "Depression; Anxiety; Obsessive-compulsive disorder",
    ),
    (
        "P08588",
        "ADRB1",
        "Beta-1 adrenergic receptor",
        "Homo sapiens",
        "GPCR",
        "G-protein coupled receptor mediating cardiac sympathetic stimulation; "
        "the target of beta-blockers.",
        "Hypertension; Angina; Heart failure; Arrhythmia",
    ),
    (
        "P12821",
        "ACE",
        "Angiotensin-converting enzyme",
        "Homo sapiens",
        "Enzyme (Protease)",
        "Zinc metallopeptidase converting angiotensin I to the vasoconstrictor "
        "angiotensin II; the target of ACE inhibitors.",
        "Hypertension; Heart failure; Diabetic nephropathy",
    ),
    (
        "P00519",
        "ABL1",
        "Tyrosine-protein kinase ABL1",
        "Homo sapiens",
        "Kinase",
        "Non-receptor tyrosine kinase; its BCR-ABL fusion drives chronic myeloid "
        "leukaemia and is the target of imatinib-class drugs.",
        "Chronic myeloid leukaemia; Acute lymphoblastic leukaemia",
    ),
    (
        "P00533",
        "EGFR",
        "Epidermal growth factor receptor",
        "Homo sapiens",
        "Kinase",
        "Receptor tyrosine kinase driving proliferation; mutated/over-expressed "
        "in many carcinomas and the target of EGFR-TKIs.",
        "Non-small-cell lung cancer; Glioblastoma; Colorectal cancer",
    ),
    (
        "P03372",
        "ESR1",
        "Estrogen receptor alpha",
        "Homo sapiens",
        "Nuclear receptor",
        "Ligand-activated transcription factor; the target of SERMs in hormone "
        "receptor-positive breast cancer.",
        "Breast cancer; Osteoporosis",
    ),
    (
        "O76074",
        "PDE5A",
        "cGMP-specific phosphodiesterase 5A",
        "Homo sapiens",
        "Enzyme (Phosphodiesterase)",
        "Hydrolyses cGMP in vascular smooth muscle; the target of erectile "
        "dysfunction and pulmonary hypertension drugs.",
        "Erectile dysfunction; Pulmonary arterial hypertension",
    ),
    (
        "P35367",
        "HRH1",
        "Histamine H1 receptor",
        "Homo sapiens",
        "GPCR",
        "G-protein coupled receptor mediating allergic responses; the target of "
        "antihistamines.",
        "Allergy; Allergic rhinitis; Urticaria",
    ),
]

# (compound name, target gene, activity_type, value_nM, units, pchembl)
ACTIVITIES: list[tuple[str, str, str, float, str, float]] = [
    # NSAIDs → COX-1 / COX-2
    ("Ibuprofen", "PTGS1", "IC50", 13000, "nM", 4.89),
    ("Ibuprofen", "PTGS2", "IC50", 7700, "nM", 5.11),
    ("Naproxen", "PTGS1", "IC50", 9000, "nM", 5.05),
    ("Naproxen", "PTGS2", "IC50", 5000, "nM", 5.30),
    ("Ketoprofen", "PTGS1", "IC50", 4000, "nM", 5.40),
    ("Ketoprofen", "PTGS2", "IC50", 9000, "nM", 5.05),
    ("Flurbiprofen", "PTGS1", "IC50", 100, "nM", 7.00),
    ("Flurbiprofen", "PTGS2", "IC50", 500, "nM", 6.30),
    ("Diclofenac", "PTGS1", "IC50", 150, "nM", 6.82),
    ("Diclofenac", "PTGS2", "IC50", 40, "nM", 7.40),
    ("Indomethacin", "PTGS1", "IC50", 100, "nM", 7.00),
    ("Indomethacin", "PTGS2", "IC50", 1000, "nM", 6.00),
    ("Aspirin", "PTGS1", "IC50", 10000, "nM", 5.00),
    ("Aspirin", "PTGS2", "IC50", 20000, "nM", 4.70),
    ("Celecoxib", "PTGS2", "IC50", 40, "nM", 7.40),
    ("Celecoxib", "PTGS1", "IC50", 15000, "nM", 4.82),
    ("Mefenamic acid", "PTGS1", "IC50", 25000, "nM", 4.60),
    ("Mefenamic acid", "PTGS2", "IC50", 1000, "nM", 6.00),
    # Statins → HMGCR
    ("Atorvastatin", "HMGCR", "IC50", 8, "nM", 8.10),
    ("Simvastatin", "HMGCR", "IC50", 11, "nM", 7.96),
    ("Pravastatin", "HMGCR", "IC50", 44, "nM", 7.36),
    ("Rosuvastatin", "HMGCR", "IC50", 5, "nM", 8.30),
    ("Fluvastatin", "HMGCR", "IC50", 28, "nM", 7.55),
    # SSRIs → SERT
    ("Fluoxetine", "SLC6A4", "Ki", 0.8, "nM", 9.10),
    ("Paroxetine", "SLC6A4", "Ki", 0.13, "nM", 9.89),
    ("Sertraline", "SLC6A4", "Ki", 0.29, "nM", 9.54),
    ("Citalopram", "SLC6A4", "Ki", 1.4, "nM", 8.85),
    ("Fluvoxamine", "SLC6A4", "Ki", 2.2, "nM", 8.66),
    # Beta-blockers → ADRB1
    ("Propranolol", "ADRB1", "Ki", 1.8, "nM", 8.74),
    ("Atenolol", "ADRB1", "Ki", 480, "nM", 6.32),
    ("Metoprolol", "ADRB1", "Ki", 250, "nM", 6.60),
    ("Bisoprolol", "ADRB1", "Ki", 26, "nM", 7.58),
    # ACE inhibitors → ACE
    ("Lisinopril", "ACE", "IC50", 1.2, "nM", 8.92),
    ("Captopril", "ACE", "IC50", 6, "nM", 8.22),
    ("Enalapril", "ACE", "IC50", 18, "nM", 7.74),
    # Kinase inhibitors
    ("Imatinib", "ABL1", "Kd", 20, "nM", 7.70),
    ("Dasatinib", "ABL1", "IC50", 0.6, "nM", 9.22),
    ("Gefitinib", "EGFR", "IC50", 2, "nM", 8.70),
    ("Erlotinib", "EGFR", "IC50", 2, "nM", 8.70),
    # SERMs → ESR1
    ("Tamoxifen", "ESR1", "Ki", 10, "nM", 8.00),
    ("Raloxifene", "ESR1", "IC50", 7, "nM", 8.15),
    # PDE5 inhibitors
    ("Sildenafil", "PDE5A", "IC50", 3.5, "nM", 8.46),
    ("Tadalafil", "PDE5A", "IC50", 5, "nM", 8.30),
    # Antihistamines → HRH1
    ("Diphenhydramine", "HRH1", "Ki", 16, "nM", 7.80),
    ("Loratadine", "HRH1", "Ki", 120, "nM", 6.92),
    ("Cetirizine", "HRH1", "Ki", 6, "nM", 8.22),
]


def get_demo_dataset() -> dict:
    """Return the demo dataset as structured lists."""
    return {"compounds": COMPOUNDS, "targets": TARGETS, "activities": ACTIVITIES}
