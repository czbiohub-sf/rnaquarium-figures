"""
Zebrafish tissue-specific marker genes reference dictionary.

Curated from published literature, ZFIN database, and validated transgenic
reporter lines. Each tissue maps to a list of dicts with fields:
    gene      - zebrafish gene symbol (ZFIN nomenclature)
    alias     - common alternative names or human ortholog
    role      - brief description of the marker's use
    transgenic- representative transgenic line(s) using this promoter
    refs      - PubMed IDs or DOIs of key validation papers

Sources (reviews that informed the compilation):
  - Choi et al. 2021, Lab Anim Res 37:28       (PMC8424172)
  - Gut et al. 2017, Development 144:4116
  - ZFIN (https://zfin.org)
  - Addgene Zebrafish Collection
"""

ZEBRAFISH_TISSUE_MARKERS = {

    # =========================================================================
    # CARDIOVASCULAR SYSTEM
    # =========================================================================
    "Cardiovascular System": [
        # --- Heart / Cardiomyocytes ---
        {
            "gene": "myl7",
            "alias": "cmlc2 (cardiac myosin light chain 2)",
            "role": "Pan-cardiac marker; most widely used heart-specific promoter in zebrafish",
            "transgenic": "Tg(myl7:EGFP), Tg(myl7:DsRed)",
            "refs": "PMID:12867029, PMID:12667537",
        },
        {
            "gene": "nkx2.5",
            "alias": "tinman ortholog",
            "role": "Cardiac progenitor transcription factor; marks first and second heart field",
            "transgenic": "Tg(nkx2.5:ZsYellow), Tg(nkx2.5:EGFP)",
            "refs": "PMID:23426627, PMID:29415608",
        },
        {
            "gene": "nppa",
            "alias": "ANF / atrial natriuretic peptide A",
            "role": "Atrial / AVC marker; distinguishes chamber identity",
            "transgenic": "Tg(nppa:mCherry)",
            "refs": "PMID:29784673",
        },
        {
            "gene": "vmhc",
            "alias": "ventricular myosin heavy chain",
            "role": "Ventricular cardiomyocyte-specific marker",
            "transgenic": "ISH probe, used with myl7 double labeling",
            "refs": "PMID:10686600",
        },
        {
            "gene": "amhc",
            "alias": "atrial myosin heavy chain (myh6)",
            "role": "Atrial cardiomyocyte-specific marker",
            "transgenic": "ISH probe",
            "refs": "PMID:11734856",
        },
        {
            "gene": "tnnt2a",
            "alias": "cardiac troponin T type 2a (silent heart)",
            "role": "Cardiac muscle contraction; mutations cause heart arrest/dilated cardiomyopathy",
            "transgenic": "Mutant: sih; CRISPR models",
            "refs": "PMID:11734856, PMID:31813949",
        },
        # --- Vascular endothelium ---
        {
            "gene": "kdrl",
            "alias": "flk1 / VEGFR2 (kdr-like)",
            "role": "Pan-endothelial marker; labels all vascular endothelial cells",
            "transgenic": "Tg(kdrl:EGFP), Tg(kdrl:mCherry), Tg(kdrl:Cre)",
            "refs": "PMID:16452098, PMID:17537913",
        },
        {
            "gene": "fli1a",
            "alias": "friend leukemia integration 1a",
            "role": "Endothelial and skeletogenic precursor marker; very widely used",
            "transgenic": "Tg(fli1a:EGFP)y1, Tg(fli1a:nEGFP)",
            "refs": "PMID:11053088, PMID:12867029",
        },
        {
            "gene": "tie2",
            "alias": "tek / angiopoietin receptor",
            "role": "Vascular endothelial cells",
            "transgenic": "Tg(Tie2:EGFP)",
            "refs": "PMID:17537913",
        },
        {
            "gene": "gata4",
            "alias": "GATA binding protein 4",
            "role": "Cardiac transcription factor; marks cardiomyocyte progenitors and epicardium",
            "transgenic": "Tg(gata4:EGFP)",
            "refs": "PMID:20707672",
        },
        # --- Smooth muscle / mural cells ---
        {
            "gene": "acta2",
            "alias": "alpha smooth muscle actin (alphaSMA)",
            "role": "Earliest mural cell marker; labels vascular smooth muscle and visceral smooth muscle",
            "transgenic": "Tg(acta2:EGFP), Tg(acta2:mCherry)",
            "refs": "PMID:24594685",
        },
        {
            "gene": "tagln",
            "alias": "transgelin / SM22alpha",
            "role": "Early smooth muscle marker; co-expressed with acta2",
            "transgenic": "Co-localizes with Tg(acta2:EGFP)",
            "refs": "PMID:24594685",
        },
        {
            "gene": "myh11",
            "alias": "smooth muscle myosin heavy chain",
            "role": "Mature smooth muscle marker; first detected in intestinal smooth muscle ~60 hpf",
            "transgenic": "ISH / immunostaining",
            "refs": "PMID:24397376",
        },
    ],

    # =========================================================================
    # NERVOUS SYSTEM
    # =========================================================================
    "Nervous System": [
        # --- Pan-neuronal ---
        {
            "gene": "elavl3",
            "alias": "HuC",
            "role": "Pan-neuronal marker; labels most post-mitotic neurons (motor and sensory)",
            "transgenic": "Tg(elavl3:EGFP), Tg(elavl3:Kaede), Tg(elavl3:GCaMP3)",
            "refs": "PMID:11053088, PMID:21951526",
        },
        {
            "gene": "elavl4",
            "alias": "HuD",
            "role": "Post-mitotic neuronal marker, partially overlapping with elavl3",
            "transgenic": "ISH probe",
            "refs": "PMID:9142858",
        },
        {
            "gene": "neurod1",
            "alias": "NeuroD / neurogenic differentiation factor 1",
            "role": "Early neural commitment and differentiation; proneural transcription factor",
            "transgenic": "Tg(neurod1:EGFP), neurod1:Cre CRISPR knock-in",
            "refs": "PMID:33547356",
        },
        {
            "gene": "tubb5",
            "alias": "beta-tubulin / neural beta-tubulin",
            "role": "Pan-neuronal cytoskeletal marker",
            "transgenic": "ISH / immunostaining (anti-acetylated tubulin)",
            "refs": "PMID:8625821",
        },
        {
            "gene": "gap43",
            "alias": "growth-associated protein 43",
            "role": "Axon growth marker; labels regenerating and growing neurons",
            "transgenic": "ISH / antibody",
            "refs": "PMID:10521381",
        },
        # --- Glial cells ---
        {
            "gene": "gfap",
            "alias": "glial fibrillary acidic protein",
            "role": "Radial glia / astroglia marker in brain, spinal cord, and retina; also marks neural stem cells",
            "transgenic": "Tg(gfap:GFP), Tg(gfap:NTR-mCherry), Tg(gfap:Cre-ERT2)",
            "refs": "PMID:16765104, PMID:19161226",
        },
        {
            "gene": "mbpa",
            "alias": "mbp (myelin basic protein a)",
            "role": "Oligodendrocyte / myelinating cell marker; marks myelin sheaths and Schwann cells",
            "transgenic": "Tg(mbp:EGFP), Tg(mbp:EGFP-CAAX), Tg(mbpa:GAL4-VP16)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "olig2",
            "alias": "oligodendrocyte lineage transcription factor 2",
            "role": "Motor neuron progenitors, oligodendrocyte lineage, radial glia of pMN domain",
            "transgenic": "Tg(olig2:EGFP) (BAC transgenic)",
            "refs": "PMID:15068796",
        },
        {
            "gene": "sox10",
            "alias": "SRY-box transcription factor 10",
            "role": "Neural crest-derived cells (Schwann cells, DRG, oligodendrocyte lineage)",
            "transgenic": "Tg(-4.9sox10:EGFP), Tg(sox10:GAL4-VP16)",
            "refs": "PMID:12586706",
        },
        {
            "gene": "nestin",
            "alias": "nes (neuroepithelial stem cell marker)",
            "role": "Neural stem cells; expressed in proliferating progenitors",
            "transgenic": "Tg(-3.9nestin:GFP), TgBAC(nes:EGFP)",
            "refs": "PMID:19161226",
        },
        # --- Specific neuronal subtypes ---
        {
            "gene": "isl1",
            "alias": "islet1 (isl1a)",
            "role": "Cranial motor neurons; early marker for motor neuron identity",
            "transgenic": "Tg(isl1:GFP) (Isl1-GFP)",
            "refs": "PMID:11053088",
        },
        {
            "gene": "slc6a3",
            "alias": "dat (dopamine transporter)",
            "role": "Dopaminergic neurons",
            "transgenic": "Tg(slc6a3:EGFP), Tg(slc6a3:CFP-NTR)",
            "refs": "PMID:22615375",
        },
        {
            "gene": "vglut2a",
            "alias": "slc17a6b (vesicular glutamate transporter)",
            "role": "Glutamatergic neurons",
            "transgenic": "Tg(vglut2a:loxP-DsRed-loxP-GFP)",
            "refs": "PMID:19553430",
        },
        {
            "gene": "gad1b",
            "alias": "glutamic acid decarboxylase 67 (GAD67)",
            "role": "GABAergic neurons",
            "transgenic": "Tg(gad1b:GFP)",
            "refs": "PMID:24694998",
        },
        {
            "gene": "phox2b",
            "alias": "paired-like homeobox 2b",
            "role": "Enteric nervous system and autonomic neurons",
            "transgenic": "Tg(-8.3bphox2b:Kaede)",
            "refs": "PMID:17960624",
        },
    ],

    # =========================================================================
    # HEMATOPOIETIC SYSTEM
    # =========================================================================
    "Hematopoietic System": [
        # --- Erythroid ---
        {
            "gene": "gata1a",
            "alias": "GATA1 (erythroid transcription factor)",
            "role": "Erythroid lineage-specific marker; earliest marker of red blood cell commitment",
            "transgenic": "Tg(gata1a:GFP), Tg(gata1a:DsRed)",
            "refs": "PMID:10440861, PMID:15509524",
        },
        {
            "gene": "hbbe1.1",
            "alias": "hemoglobin beta embryonic 1.1",
            "role": "Embryonic erythrocyte marker; expressed in primitive erythrocytes",
            "transgenic": "ISH probe",
            "refs": "PMID:9142858",
        },
        {
            "gene": "hbae1.1",
            "alias": "hemoglobin alpha embryonic 1.1",
            "role": "Embryonic erythrocyte marker",
            "transgenic": "ISH probe",
            "refs": "ZFIN:ZDB-GENE-980526-80",
        },
        # --- Myeloid / neutrophils ---
        {
            "gene": "mpx",
            "alias": "mpo (myeloperoxidase)",
            "role": "Neutrophil-specific marker",
            "transgenic": "Tg(mpx:GFP)",
            "refs": "PMID:16651656",
        },
        {
            "gene": "lyz",
            "alias": "lysozyme C",
            "role": "Myeloid-specific marker (neutrophils and monocytes)",
            "transgenic": "Tg(lyz:EGFP), Tg(lyz:DsRed)",
            "refs": "PMID:17510218",
        },
        # --- Macrophages ---
        {
            "gene": "mpeg1.1",
            "alias": "mpeg1 (macrophage expressed gene 1)",
            "role": "Macrophage-specific marker; nonoverlapping with neutrophil markers mpx/lyz",
            "transgenic": "Tg(mpeg1:EGFP), Tg(mpeg1:mCherry)",
            "refs": "PMID:21297037",
        },
        {
            "gene": "mfap4",
            "alias": "microfibril associated protein 4",
            "role": "Macrophage marker; alternative to mpeg1",
            "transgenic": "Tg(mfap4:tdTomato)",
            "refs": "PMID:25510288",
        },
        # --- Lymphoid ---
        {
            "gene": "rag1",
            "alias": "recombination activating gene 1",
            "role": "Lymphocyte marker (T cells and immature B cells in thymus)",
            "transgenic": "Tg(rag1:GFP)",
            "refs": "PMID:10371508",
        },
        {
            "gene": "rag2",
            "alias": "recombination activating gene 2",
            "role": "Lymphocyte marker; labels immature lymphocytes",
            "transgenic": "Tg(rag2:GFP)",
            "refs": "PMID:15040818",
        },
        {
            "gene": "lck",
            "alias": "lymphocyte-specific protein tyrosine kinase",
            "role": "T cell-specific marker",
            "transgenic": "Tg(lck:GFP)",
            "refs": "PMID:14963217",
        },
        # --- Thrombocytes ---
        {
            "gene": "itga2b",
            "alias": "cd41 (integrin alpha 2b)",
            "role": "Thrombocyte marker; also labels early HSPCs",
            "transgenic": "Tg(-6.0itga2b:EGFP) / Tg(cd41:GFP)",
            "refs": "PMID:15750594, PMID:19755505",
        },
        # --- HSPCs (hematopoietic stem/progenitor cells) ---
        {
            "gene": "runx1",
            "alias": "Runx1 / AML1",
            "role": "Earliest marker of definitive HSPCs; required for hemogenic endothelium",
            "transgenic": "Tg(Mmu.Runx1:EGFP), Tg(Mmu.Runx1:NLS-mCherry)",
            "refs": "PMID:19652178",
        },
        {
            "gene": "cmyb",
            "alias": "myb (c-Myb proto-oncogene)",
            "role": "HSC marker; downstream of runx1 in definitive hematopoiesis",
            "transgenic": "TgPAC(myb:2xmyb-EGFP)",
            "refs": "PMID:17473174",
        },
        {
            "gene": "lmo2",
            "alias": "LIM domain only 2",
            "role": "Hemangioblast / HSC marker; marks hematopoietic and endothelial progenitors",
            "transgenic": "Tg(lmo2:DsRed)",
            "refs": "PMID:12842908",
        },
        {
            "gene": "gata2b",
            "alias": "GATA binding protein 2b",
            "role": "Hemogenic endothelium-specific; precursor of definitive HSCs",
            "transgenic": "TgBAC(gata2b:KALTA4)",
            "refs": "PMID:25788672",
        },
    ],

    # =========================================================================
    # MUSCULAR SYSTEM
    # =========================================================================
    "Muscular System": [
        # --- Fast skeletal muscle ---
        {
            "gene": "mylz2",
            "alias": "mylpfa (myosin light chain, phosphorylatable, fast skeletal muscle a)",
            "role": "Fast skeletal muscle-specific marker; widely used promoter for skeletal muscle",
            "transgenic": "Tg(mylz2:GFP)gz8",
            "refs": "PMID:12701095",
        },
        {
            "gene": "acta1a",
            "alias": "skeletal alpha-actin (acta1)",
            "role": "Skeletal muscle actin; labels myofibers (myofibrils and sarcolemma)",
            "transgenic": "Tg(acta1:lifeact-GFP), Tg(acta1:mCherryCAAX)",
            "refs": "PMID:28761223",
        },
        {
            "gene": "mylpfa",
            "alias": "myosin light polypeptide fast skeletal a",
            "role": "Fast muscle fiber-specific marker",
            "transgenic": "Tg(mylpfa:mCherry)",
            "refs": "PMID:12701095",
        },
        # --- Slow skeletal muscle ---
        {
            "gene": "smyhc1",
            "alias": "slow myosin heavy chain 1 (myh7-like)",
            "role": "Slow skeletal muscle-specific marker; slow-twitch fiber identity",
            "transgenic": "Tg(smyhc1:GFP)",
            "refs": "PMID:19883647",
        },
        {
            "gene": "tnnt1",
            "alias": "troponin T1 (slow skeletal muscle isoform)",
            "role": "Slow muscle fiber marker",
            "transgenic": "ISH / immunostaining (F59 antibody)",
            "refs": "PMID:14568581",
        },
        # --- Muscle progenitors / satellite cells ---
        {
            "gene": "myod1",
            "alias": "myoD (myogenic differentiation 1)",
            "role": "Myogenic transcription factor; marks skeletal muscle progenitors",
            "transgenic": "Tg(myod1:GFP)",
            "refs": "PMID:12524040",
        },
        {
            "gene": "myf5",
            "alias": "myogenic factor 5",
            "role": "Early myogenic progenitor marker; precedes myod activation",
            "transgenic": "Tg(myf5:GFP)",
            "refs": "PMID:12524040",
        },
        {
            "gene": "myog",
            "alias": "myogenin",
            "role": "Late myogenic differentiation factor; marks differentiating muscle",
            "transgenic": "Tg(myog:GFP)",
            "refs": "PMID:12524040",
        },
        {
            "gene": "pax7a",
            "alias": "paired box 7a",
            "role": "Muscle satellite/stem cell marker; quiescent and activated progenitors",
            "transgenic": "Tg(pax7a:GFP)",
            "refs": "PMID:22406451",
        },
        {
            "gene": "desma",
            "alias": "desmin a",
            "role": "Intermediate filament in muscle; pan-muscle structural protein",
            "transgenic": "Tg(desma:EGFP)",
            "refs": "PMID:28085148",
        },
    ],

    # =========================================================================
    # LIVER AND BILIARY SYSTEM
    # =========================================================================
    "Liver and Biliary System": [
        # --- Hepatocytes ---
        {
            "gene": "fabp10a",
            "alias": "lfabp (liver fatty acid binding protein 10a)",
            "role": "Hepatocyte-specific marker; most widely used liver promoter in zebrafish",
            "transgenic": "Tg(-2.8fabp10a:EGFP), Tg(fabp10a:DsRed)",
            "refs": "PMID:15694379, PMID:21953180",
        },
        {
            "gene": "hhex",
            "alias": "hematopoietically expressed homeobox",
            "role": "Hepatoblast specification marker; earliest liver progenitor marker at ~22 hpf",
            "transgenic": "ISH probe",
            "refs": "PMID:11597192",
        },
        {
            "gene": "prox1a",
            "alias": "prospero homeobox 1a",
            "role": "Hepatoblast specification and hepatocyte differentiation marker",
            "transgenic": "ISH probe; Tg(prox1a:GFP) labels liver and lymphatics",
            "refs": "PMID:11832215",
        },
        {
            "gene": "cp",
            "alias": "ceruloplasmin",
            "role": "Mature hepatocyte marker; expressed during hepatocyte differentiation (30-48 hpf)",
            "transgenic": "ISH probe",
            "refs": "PMID:18585888",
        },
        {
            "gene": "gc",
            "alias": "group-specific component / vitamin D binding protein",
            "role": "Hepatocyte differentiation marker",
            "transgenic": "ISH probe",
            "refs": "PMID:18585888",
        },
        {
            "gene": "bhmt",
            "alias": "betaine-homocysteine S-methyltransferase",
            "role": "Hepatocyte-specific metabolic marker",
            "transgenic": "ISH probe",
            "refs": "PMID:18787121",
        },
        # --- Biliary / cholangiocytes ---
        {
            "gene": "krt18a.1",
            "alias": "keratin 18 / cytokeratin 18",
            "role": "Intrahepatic and extrahepatic biliary epithelial cell marker",
            "transgenic": "Tg(krt18:GFP)",
            "refs": "PMID:19900452",
        },
        {
            "gene": "tp1",
            "alias": "Notch reporter (EPV.TP1-Mmu.Hbb)",
            "role": "Intrahepatic biliary cell marker (Notch-responsive)",
            "transgenic": "Tg(EPV.TP1-Mmu.Hbb:EGFP)",
            "refs": "PMID:17070912",
        },
        {
            "gene": "anxa4",
            "alias": "annexin A4",
            "role": "Biliary epithelial cell marker",
            "transgenic": "Tg(anxa4:GFP)",
            "refs": "PMID:24436323",
        },
        # --- Hepatic stellate cells ---
        {
            "gene": "hand2",
            "alias": "heart and neural crest derivatives expressed 2",
            "role": "Hepatic stellate cell marker",
            "transgenic": "TgBAC(hand2:EGFP)",
            "refs": "PMID:21953180",
        },
    ],

    # =========================================================================
    # DIGESTIVE SYSTEM
    # =========================================================================
    "Digestive System": [
        # --- Intestine ---
        {
            "gene": "fabp2",
            "alias": "I-FABP (intestinal fatty acid binding protein)",
            "role": "Intestinal epithelium marker; 4.5kb promoter drives gut-specific expression",
            "transgenic": "Tg(fabp2:RFP), Tg(fabp2:GFP)",
            "refs": "PMID:15232758, PMID:21953180",
        },
        {
            "gene": "vil1",
            "alias": "villin 1",
            "role": "Intestinal brush border marker; small intestine-enriched",
            "transgenic": "ISH probe",
            "refs": "PMID:19884307",
        },
        {
            "gene": "fabp6",
            "alias": "ileal fatty acid binding protein",
            "role": "Ileum-like region marker (mid-intestine)",
            "transgenic": "ISH probe",
            "refs": "PMID:19884307",
        },
        {
            "gene": "cldn15la",
            "alias": "claudin 15-like a",
            "role": "Intestinal epithelial cell marker (bulb to posterior intestine)",
            "transgenic": "TgBAC(cldn15la:GFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "apoa1a",
            "alias": "apolipoprotein A-Ia",
            "role": "Anterior intestine marker (jejunum-like)",
            "transgenic": "ISH probe",
            "refs": "PMID:19884307",
        },
        # --- Endoderm / gut specification ---
        {
            "gene": "foxa3",
            "alias": "forkhead box A3 (HNF3gamma)",
            "role": "Endoderm marker; drives expression in developing gut, liver, and pancreas",
            "transgenic": "Tg(foxa3:GFP)",
            "refs": "PMID:12890018",
        },
        {
            "gene": "sox17",
            "alias": "SRY-box 17",
            "role": "Endoderm and pharyngeal endoderm marker",
            "transgenic": "Tg(sox17:GFP)",
            "refs": "PMID:11734856",
        },
        # --- Pancreas exocrine ---
        {
            "gene": "ela3l",
            "alias": "elastase A / elaA (pancreatic elastase 3 like)",
            "role": "Pancreatic exocrine cell-specific marker",
            "transgenic": "Tg(ela3l:EGFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "try",
            "alias": "trypsin",
            "role": "Exocrine pancreas marker; digestive enzyme",
            "transgenic": "Tg(ptf1a:EGFP) labels exocrine progenitors",
            "refs": "PMID:12167406",
        },
        {
            "gene": "ptf1a",
            "alias": "pancreas transcription factor 1a",
            "role": "Exocrine pancreas progenitor transcription factor",
            "transgenic": "Tg(ptf1a:EGFP)",
            "refs": "PMID:12167406",
        },
        # --- Pancreas endocrine (see also Endocrine System) ---
        {
            "gene": "ins",
            "alias": "insulin (preproinsulin)",
            "role": "Pancreatic beta-cell marker",
            "transgenic": "Tg(ins:DsRed), Tg(ins:GFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "gcga",
            "alias": "glucagon a",
            "role": "Pancreatic alpha-cell marker",
            "transgenic": "Tg(gcga:GFP), Tg(gcga:Cre)",
            "refs": "PMID:21953180, PMID:37516089",
        },
    ],

    # =========================================================================
    # REPRODUCTIVE SYSTEM
    # =========================================================================
    "Reproductive System": [
        {
            "gene": "ddx4",
            "alias": "vasa (DEAD-box helicase 4)",
            "role": "Conserved germ cell marker; labels primordial germ cells (PGCs) and all germ lineage",
            "transgenic": "Tg(ddx4:ddx4-EGFP) / Tg(vasa:vasa-EGFP)",
            "refs": "PMID:11463859",
        },
        {
            "gene": "piwil1",
            "alias": "ziwi (piwi-like RNA-mediated gene silencing 1)",
            "role": "Germ cell marker; expressed in PGCs and gamete progenitors",
            "transgenic": "Tg(piwil1:EGFP), Tg(piwil1:EGFP-UTRnanos3)",
            "refs": "PMID:30814499",
        },
        {
            "gene": "nanos3",
            "alias": "nanos C2HC-type zinc finger 3",
            "role": "PGC specification and maintenance; 3-prime UTR used for germline-specific mRNA stability",
            "transgenic": "nanos3 3-prime-UTR used in Tg(piwil1:EGFP-UTRnanos3)",
            "refs": "PMID:30814499",
        },
        {
            "gene": "dnd1",
            "alias": "dead end 1 (RNA-binding protein)",
            "role": "Essential for PGC survival; loss causes all-male fish (no germ cells)",
            "transgenic": "dnd1 morpholino; mutant lines",
            "refs": "PMID:14973291",
        },
        {
            "gene": "gsdf",
            "alias": "gonadal somatic cell derived factor",
            "role": "Sertoli and granulosa cell marker; gonadal somatic cells",
            "transgenic": "Tg(gsdf:EGFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "amh",
            "alias": "anti-Muellerian hormone",
            "role": "Sertoli cell marker in testis; sex differentiation",
            "transgenic": "ISH probe",
            "refs": "PMID:12871295",
        },
        {
            "gene": "cyp19a1a",
            "alias": "aromatase (gonadal)",
            "role": "Ovary-specific marker; converts androgens to estrogens",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:15950195",
        },
        {
            "gene": "piwil2",
            "alias": "zili (piwi-like 2)",
            "role": "Germ cell marker; required for meiosis and transposon silencing",
            "transgenic": "Mutant lines",
            "refs": "PMID:20360737",
        },
    ],

    # =========================================================================
    # RENAL SYSTEM
    # =========================================================================
    "Renal System": [
        # --- Podocytes / Glomerulus ---
        {
            "gene": "nphs2",
            "alias": "podocin",
            "role": "Podocyte-specific marker; glomerular epithelial cells",
            "transgenic": "Tg(-2.5nphs2:EGFP), Tg(nphs2:GAL4-VP16)",
            "refs": "PMID:22440901",
        },
        {
            "gene": "wt1a",
            "alias": "Wilms tumor 1a",
            "role": "Podocyte progenitor marker; required for glomerular formation",
            "transgenic": "Tg(wt1b:EGFP) (wt1b paralog used more widely for transgenic)",
            "refs": "PMID:17437786",
        },
        {
            "gene": "wt1b",
            "alias": "Wilms tumor 1b",
            "role": "Pronephric kidney marker; glomerulus and tubule specification",
            "transgenic": "Tg(wt1b:EGFP)",
            "refs": "PMID:17437786",
        },
        {
            "gene": "nphs1",
            "alias": "nephrin",
            "role": "Podocyte slit diaphragm marker; ortholog of NPHS1",
            "transgenic": "ISH probe",
            "refs": "PMID:22440901",
        },
        # --- Tubule ---
        {
            "gene": "cdh17",
            "alias": "kidney-specific cadherin (cadherin 17)",
            "role": "Pronephric tubule and duct marker (whole length); excludes glomeruli",
            "transgenic": "Tg(cdh17:EGFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "slc20a1a",
            "alias": "sodium-phosphate co-transporter",
            "role": "Proximal tubule marker",
            "transgenic": "ISH probe",
            "refs": "PMID:20045474",
        },
        {
            "gene": "enpep",
            "alias": "glutamyl aminopeptidase",
            "role": "Pronephric tubules, ducts, and podocyte-like cells of glomeruli",
            "transgenic": "Tg(enpep:EGFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "pax2a",
            "alias": "paired box 2a",
            "role": "Kidney progenitor marker; co-expressed with wt1a/wt1b in early nephron precursors",
            "transgenic": "Tg(pax2a:GFP)",
            "refs": "PMID:10207054",
        },
        {
            "gene": "lhx1a",
            "alias": "LIM homeobox 1a",
            "role": "Pronephric mesoderm and tubule specification marker",
            "transgenic": "Tg(lhx1a:EGFP)",
            "refs": "PMID:21953180",
        },
    ],

    # =========================================================================
    # SENSORY SYSTEM (eye / ear / lateral line)
    # =========================================================================
    "Sensory System": [
        # --- Eye: retinal progenitors ---
        {
            "gene": "rx1",
            "alias": "retinal homeobox gene 1 (rax)",
            "role": "Retinal progenitor marker; defines region of forebrain giving rise to retina",
            "transgenic": "ISH probe",
            "refs": "PMID:10529420",
        },
        {
            "gene": "pax6a",
            "alias": "paired box 6a",
            "role": "Eye field specification; master regulator of eye development",
            "transgenic": "Tg(pax6a:GFP)",
            "refs": "PMID:17274609",
        },
        {
            "gene": "vsx2",
            "alias": "visual system homeobox 2 (Chx10)",
            "role": "Retinal progenitor cells and bipolar cells",
            "transgenic": "Tg(vsx2:GFP)",
            "refs": "PMID:15031319",
        },
        # --- Eye: retinal ganglion cells ---
        {
            "gene": "atoh7",
            "alias": "ath5 / lakritz",
            "role": "Retinal ganglion cell (RGC) specification; required for all RGC genesis",
            "transgenic": "Tg(atoh7:GFP); mutant: lakritz (lak)",
            "refs": "PMID:11500489, PMID:32597793",
        },
        {
            "gene": "isl2b",
            "alias": "islet2b",
            "role": "Retinal ganglion cell-specific marker with robust adult expression",
            "transgenic": "Tg(isl2b:GFP)",
            "refs": "PMID:33357413",
        },
        {
            "gene": "pou4f3",
            "alias": "Brn3c",
            "role": "Hair cell marker (ear and lateral line); also labels RGC subset",
            "transgenic": "Tg(pou4f3:gap43-GFP)",
            "refs": "PMID:21482785",
        },
        # --- Eye: photoreceptors ---
        {
            "gene": "crx",
            "alias": "cone-rod homeobox",
            "role": "Photoreceptor transcription factor; marks all photoreceptors",
            "transgenic": "Tg(crx:mYFP)",
            "refs": "PMID:25999792",
        },
        {
            "gene": "rho",
            "alias": "rhodopsin (rod opsin)",
            "role": "Rod photoreceptor marker",
            "transgenic": "Tg(rho:EGFP)",
            "refs": "PMID:14561723",
        },
        {
            "gene": "opn1sw1",
            "alias": "UV cone opsin (SWS1)",
            "role": "UV-sensitive cone photoreceptor marker",
            "transgenic": "Tg(opn1sw1:GFP)",
            "refs": "PMID:26260523",
        },
        {
            "gene": "opn1sw2",
            "alias": "blue cone opsin (SWS2)",
            "role": "Blue-sensitive cone photoreceptor marker",
            "transgenic": "ISH probe",
            "refs": "PMID:26260523",
        },
        {
            "gene": "opn1mw1",
            "alias": "green cone opsin (RH2)",
            "role": "Green-sensitive cone photoreceptor marker",
            "transgenic": "ISH probe",
            "refs": "PMID:26260523",
        },
        {
            "gene": "opn1lw1",
            "alias": "red cone opsin (LWS)",
            "role": "Red-sensitive cone photoreceptor marker",
            "transgenic": "ISH probe",
            "refs": "PMID:26260523",
        },
        # --- Eye: lens ---
        {
            "gene": "cryaa",
            "alias": "alpha-crystallin A",
            "role": "Lens-specific marker; unique lens expression in larval and adult fish",
            "transgenic": "Tg(cryaa:EGFP), Tg(cryaa:Cre)",
            "refs": "PMID:38534198",
        },
        # --- Ear / inner ear ---
        {
            "gene": "atoh1a",
            "alias": "atonal bHLH transcription factor 1a",
            "role": "Hair cell specification in inner ear and lateral line neuromasts",
            "transgenic": "Tg(atoh1a:GFP)",
            "refs": "PMID:15504907",
        },
        {
            "gene": "myo6b",
            "alias": "myosin VI b",
            "role": "Hair cell marker (inner ear and lateral line); cytoskeletal component",
            "transgenic": "ISH / scRNA-seq marker",
            "refs": "PMID:36583561",
        },
        {
            "gene": "eya1",
            "alias": "eyes absent homolog 1",
            "role": "Otic vesicle and lateral line marker; part of Pax-Six-Eya-Dach network",
            "transgenic": "ISH probe; mutant: dog-eared",
            "refs": "PMID:17035528",
        },
        {
            "gene": "six1b",
            "alias": "sine oculis homeobox 1b",
            "role": "Preplacodal ectoderm and otic/lateral line placode marker",
            "transgenic": "ISH probe",
            "refs": "PMID:17035528",
        },
        # --- Lateral line ---
        {
            "gene": "cldnb",
            "alias": "claudin b",
            "role": "Neuromast and lateral line primordium marker (whole neuromast + interneuromast cells)",
            "transgenic": "Tg(cldnb:lynGFP)",
            "refs": "PMID:16049015",
        },
        {
            "gene": "cxcr4b",
            "alias": "C-X-C motif chemokine receptor 4b",
            "role": "Lateral line primordium leading zone chemokine receptor; migration guidance",
            "transgenic": "Tg(cxcr4b:mRFP)",
            "refs": "PMID:15198164",
        },
    ],

    # =========================================================================
    # ENDOCRINE SYSTEM
    # =========================================================================
    "Endocrine System": [
        # --- Pancreatic endocrine ---
        {
            "gene": "ins",
            "alias": "insulin (beta cells)",
            "role": "Pancreatic beta-cell marker; most used endocrine pancreas promoter",
            "transgenic": "Tg(ins:DsRed), Tg(ins:GFP), Tg(ins:NTR-mCherry)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "gcga",
            "alias": "glucagon a (alpha cells)",
            "role": "Pancreatic alpha-cell marker",
            "transgenic": "Tg(gcga:GFP), Tg(gcga:Cre)",
            "refs": "PMID:37516089",
        },
        {
            "gene": "sst1.1",
            "alias": "somatostatin 1.1 (delta cells)",
            "role": "Pancreatic delta-cell marker; also expressed in brain neuroendocrine cells",
            "transgenic": "Tg(sst1.1:GFP)",
            "refs": "PMID:32661042",
        },
        {
            "gene": "ghrl",
            "alias": "ghrelin (epsilon cells)",
            "role": "Ghrelin-producing epsilon-cell marker in pancreas",
            "transgenic": "ISH / scRNA-seq marker",
            "refs": "PMID:32661042",
        },
        {
            "gene": "mnx1",
            "alias": "motor neuron and pancreas homeobox 1 (Hb9)",
            "role": "Beta-cell marker; transcription factor for beta-cell fate",
            "transgenic": "Tg(mnx1:GFP)",
            "refs": "PMID:21989909",
        },
        {
            "gene": "pdx1",
            "alias": "pancreatic and duodenal homeobox 1",
            "role": "Pancreatic progenitor and beta-cell transcription factor",
            "transgenic": "Tg(pdx1:GFP)",
            "refs": "PMID:21989909",
        },
        {
            "gene": "neurod1",
            "alias": "neurogenic differentiation 1",
            "role": "Endocrine pancreas differentiation; controls alpha/beta cell fate balance",
            "transgenic": "Tg(neurod1:EGFP)",
            "refs": "PMID:25989474",
        },
        # --- Pituitary ---
        {
            "gene": "pomca",
            "alias": "proopiomelanocortin a",
            "role": "Pituitary corticotrope and melanotrope marker",
            "transgenic": "Tg(pomc:EGFP)",
            "refs": "PMID:32434973",
        },
        {
            "gene": "gh1",
            "alias": "growth hormone 1",
            "role": "Pituitary somatotrope marker",
            "transgenic": "ISH probe",
            "refs": "PMID:35149677",
        },
        {
            "gene": "prl",
            "alias": "prolactin",
            "role": "Pituitary lactotrope marker",
            "transgenic": "ISH probe",
            "refs": "PMID:35149677",
        },
        {
            "gene": "tshba",
            "alias": "thyroid stimulating hormone beta a",
            "role": "Pituitary thyrotrope marker",
            "transgenic": "Tg(tshba:GFP) / Tg(gtshb:GFP)",
            "refs": "PMID:35149677",
        },
        # --- Thyroid ---
        {
            "gene": "tg",
            "alias": "thyroglobulin",
            "role": "Thyroid follicular cell marker",
            "transgenic": "ISH probe",
            "refs": "PMID:11756311",
        },
        {
            "gene": "nkx2.1a",
            "alias": "thyroid transcription factor 1 (TTF1)",
            "role": "Thyroid primordium specification marker",
            "transgenic": "ISH probe",
            "refs": "PMID:11756311",
        },
        # --- Interrenal (adrenal analog) ---
        {
            "gene": "cyp11a1",
            "alias": "cytochrome P450 side chain cleavage (P450scc)",
            "role": "Interrenal (adrenal cortex analog) steroidogenic cell marker",
            "transgenic": "ISH probe",
            "refs": "PMID:12842912",
        },
        {
            "gene": "star",
            "alias": "steroidogenic acute regulatory protein",
            "role": "Steroidogenic cell marker (interrenal gland)",
            "transgenic": "ISH probe",
            "refs": "PMID:12842912",
        },
    ],

    # =========================================================================
    # SKELETAL ELEMENT
    # =========================================================================
    "Skeletal Element": [
        # --- Osteoblasts / Bone ---
        {
            "gene": "sp7",
            "alias": "osterix (Osx)",
            "role": "Osteoblast-specific transcription factor; essential for bone mineralization",
            "transgenic": "Tg(Ola.Sp7:nlsGFP), Tg(sp7:EGFP)",
            "refs": "PMID:21953180, PMID:31574273",
        },
        {
            "gene": "bglap",
            "alias": "osteocalcin (bone gamma-carboxyglutamate protein)",
            "role": "Late / mature osteoblast marker; marks mineralized bone",
            "transgenic": "Tg(Ola.Bglap:EGFP)",
            "refs": "PMID:21953180, PMID:19288476",
        },
        {
            "gene": "runx2a",
            "alias": "runt-related transcription factor 2a",
            "role": "Early osteoblast transcription factor; master regulator of osteoblast differentiation",
            "transgenic": "Tg(Hsa.RUNX2-Mmu.Fos:EGFP)",
            "refs": "PMID:21953180, PMID:19288476",
        },
        {
            "gene": "runx2b",
            "alias": "runt-related transcription factor 2b",
            "role": "Osteoblast differentiation; duplicate of runx2a with subfunctionalization",
            "transgenic": "ISH / mutant studies",
            "refs": "PMID:36103880",
        },
        {
            "gene": "alpl",
            "alias": "alkaline phosphatase, liver/bone/kidney",
            "role": "Osteoblast activity marker; enzyme marker for bone formation",
            "transgenic": "Enzyme histochemistry / ISH",
            "refs": "PMID:19288476",
        },
        {
            "gene": "col1a1a",
            "alias": "collagen type I alpha 1a",
            "role": "Osteoblast matrix protein; major bone collagen",
            "transgenic": "Tg(col1a1a:EGFP)",
            "refs": "PMID:19288476",
        },
        # --- Chondrocytes / Cartilage ---
        {
            "gene": "col2a1a",
            "alias": "collagen type II alpha 1a",
            "role": "Chondrocyte marker; primary cartilage collagen",
            "transgenic": "Tg(col2a1a:EGFP), Tg(Col2a1aBAC:mCherry)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "sox9a",
            "alias": "SRY-box transcription factor 9a",
            "role": "Chondrocyte master transcription factor; required for cartilage and dermal bone",
            "transgenic": "Tg(sox9a:EGFP)",
            "refs": "PMID:21953180",
        },
        {
            "gene": "sox9b",
            "alias": "SRY-box transcription factor 9b",
            "role": "Cartilage and bone development; partially redundant with sox9a",
            "transgenic": "ISH / mutant studies",
            "refs": "PMID:15843405",
        },
        {
            "gene": "col10a1",
            "alias": "collagen type X alpha 1 (hypertrophic chondrocyte collagen)",
            "role": "Hypertrophic chondrocyte marker; marks chondrocytes undergoing endochondral ossification",
            "transgenic": "ISH probe",
            "refs": "PMID:22329103",
        },
        {
            "gene": "acan",
            "alias": "aggrecan (cartilage proteoglycan)",
            "role": "Cartilage matrix marker; major proteoglycan in cartilage",
            "transgenic": "ISH probe",
            "refs": "PMID:22329103",
        },
    ],

    # =========================================================================
    # SURFACE STRUCTURE (SKIN / EPIDERMIS)
    # =========================================================================
    "Surface Structure": [
        {
            "gene": "krt4",
            "alias": "keratin 4 (cytokeratin, type II)",
            "role": "Periderm / superficial skin marker; labels both EVL layers; most used skin promoter",
            "transgenic": "Tg(krt4:EGFP), Tg(krt4:DsRed), multiple enhancer trap lines Et(krt4:EGFP)",
            "refs": "PMID:20580760, PMID:21479254",
        },
        {
            "gene": "krt5",
            "alias": "keratin 5",
            "role": "Periderm (outer skin layer) specific at later stages; not in basal cells",
            "transgenic": "Tg(krt5:GFP)",
            "refs": "PMID:31451637",
        },
        {
            "gene": "krt8",
            "alias": "keratin 8 (cytokeratin 8)",
            "role": "Envelope layer / stratified epithelium marker (skin, mouth, pharynx, esophagus)",
            "transgenic": "Tg(krt8:EGFP)",
            "refs": "PMID:11836785",
        },
        {
            "gene": "krt18a.1",
            "alias": "keratin 18",
            "role": "Simple epithelium marker; periderm and digestive tract epithelium",
            "transgenic": "Tg(krt18:GFP)",
            "refs": "PMID:19900452",
        },
        {
            "gene": "tp63",
            "alias": "tumor protein p63 (deltaNp63)",
            "role": "Basal epidermal cell / skin stem cell marker; not in periderm/EVL",
            "transgenic": "ISH / antibody; Tg(tp63:GFP)",
            "refs": "PMID:21479254",
        },
        {
            "gene": "krtt1c19e",
            "alias": "keratin type I c19e",
            "role": "Skin epithelium marker",
            "transgenic": "ISH probe",
            "refs": "PMID:31451637",
        },
        {
            "gene": "cyt1",
            "alias": "cytokeratin type I (enveloping layer keratin)",
            "role": "EVL/periderm differentiation marker; part of Grhl/Klf17 regulatory network",
            "transgenic": "ISH probe",
            "refs": "PMID:31451637",
        },
    ],

    # =========================================================================
    # SWIM BLADDER
    # =========================================================================
    "Swim Bladder": [
        {
            "gene": "hb9",
            "alias": "mnx1 (motor neuron and pancreas homeobox 1)",
            "role": "Swim bladder epithelium marker; earliest detectable epithelial marker",
            "transgenic": "ISH probe; Tg(mnx1:GFP) also labels swim bladder epithelium",
            "refs": "PMID:19422819",
        },
        {
            "gene": "anxa5b",
            "alias": "annexin A5b",
            "role": "Swim bladder mesothelium (outer layer) marker",
            "transgenic": "ISH probe; Et(krt4:EGFP)sqet3 labels mesothelium",
            "refs": "PMID:19422819, PMID:20108353",
        },
        {
            "gene": "fgf10a",
            "alias": "fibroblast growth factor 10a",
            "role": "Swim bladder mesenchyme layer marker; analogous to lung bud FGF10",
            "transgenic": "ISH probe",
            "refs": "PMID:19422819",
        },
        {
            "gene": "pbx1a",
            "alias": "pre-B-cell leukemia homeobox 1a",
            "role": "Earliest known swim bladder marker (~28 hpf); essential for swim bladder growth",
            "transgenic": "pbx1 morpholino; ISH probe",
            "refs": "PMID:20108353",
        },
        {
            "gene": "shha",
            "alias": "sonic hedgehog signaling molecule a",
            "role": "Swim bladder epithelial signaling; marks epithelium (conserved with lung Shh)",
            "transgenic": "ISH probe",
            "refs": "PMID:19422819",
        },
        {
            "gene": "acta2",
            "alias": "alpha smooth muscle actin",
            "role": "Swim bladder smooth muscle (mesenchyme layer) marker",
            "transgenic": "Tg(acta2:EGFP) labels swim bladder mesenchyme",
            "refs": "PMID:19422819",
        },
        {
            "gene": "sftpc",
            "alias": "surfactant protein C (putative)",
            "role": "Swim bladder surfactant system; homologous to mammalian lung surfactant",
            "transgenic": "ISH probe",
            "refs": "PMID:29516002",
        },
        {
            "gene": "sox2",
            "alias": "SRY-box transcription factor 2",
            "role": "Swim bladder inflation regulation; controls swim-up behavior and inflation",
            "transgenic": "sox2 mutant studies",
            "refs": "PMID:36825027",
        },
    ],

    # =========================================================================
    # RESPIRATORY SYSTEM (gills)
    # =========================================================================
    "Respiratory System": [
        {
            "gene": "atp1b1b",
            "alias": "Na+/K+ ATPase beta 1b",
            "role": "Gill ionocyte (NaR cell) marker; ion transport in gill epithelium",
            "transgenic": "ISH / immunostaining",
            "refs": "PMID:12574125",
        },
        {
            "gene": "atp1a1a.1",
            "alias": "Na+/K+ ATPase alpha 1a.1",
            "role": "Gill ionocyte marker; Na+/K+ ATPase-rich cell identity",
            "transgenic": "ISH / immunostaining",
            "refs": "PMID:12574125",
        },
        {
            "gene": "trpv6",
            "alias": "transient receptor potential vanilloid 6 (ECaC)",
            "role": "Gill epithelial calcium channel; marks NaR ionocytes",
            "transgenic": "ISH probe",
            "refs": "PMID:15883196",
        },
        {
            "gene": "slc12a10.2",
            "alias": "NCC (Na-Cl cotransporter)",
            "role": "Gill NCC-type ionocyte marker",
            "transgenic": "ISH probe",
            "refs": "PMID:18497255",
        },
        {
            "gene": "gcm2",
            "alias": "glial cells missing transcription factor 2",
            "role": "Gill filament specification; marks developing gill arches",
            "transgenic": "ISH probe; mutant studies",
            "refs": "PMID:12533115",
        },
        {
            "gene": "foxi3a",
            "alias": "forkhead box I3a",
            "role": "Ionocyte progenitor marker; required for ionocyte specification in skin and gills",
            "transgenic": "ISH probe",
            "refs": "PMID:18339674",
        },
        {
            "gene": "hbbe1.1",
            "alias": "hemoglobin beta embryonic 1.1",
            "role": "Respiratory gas transport; embryonic hemoglobin in circulating erythrocytes",
            "transgenic": "ISH probe",
            "refs": "ZFIN:ZDB-GENE-980526-80",
        },
        {
            "gene": "hbae1.1",
            "alias": "hemoglobin alpha embryonic 1.1",
            "role": "Respiratory gas transport; embryonic hemoglobin",
            "transgenic": "ISH probe",
            "refs": "ZFIN:ZDB-GENE-980526-80",
        },
    ],

    # =========================================================================
    # ADIPOSE TISSUE
    # =========================================================================
    "Adipose Tissue": [
        {
            "gene": "fabp11a",
            "alias": "fatty acid binding protein 11a (adipocyte FABP)",
            "role": "Adipocyte-specific marker; 1.5kb promoter sufficient for adipocyte-specific expression",
            "transgenic": "Tg(fabp11a:EGFP)",
            "refs": "PMID:33826568",
        },
        {
            "gene": "pparg",
            "alias": "peroxisome proliferator-activated receptor gamma",
            "role": "Master regulator of adipogenesis; colocalizes with visceral lipid droplets",
            "transgenic": "ISH probe; co-localizes with fabp11a",
            "refs": "PMID:21951526, PMID:19664743",
        },
        {
            "gene": "cebpa",
            "alias": "CCAAT/enhancer binding protein alpha",
            "role": "Adipocyte differentiation transcription factor",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:19664743",
        },
        {
            "gene": "lpl",
            "alias": "lipoprotein lipase",
            "role": "Adipocyte lipolysis marker; lipid metabolism enzyme in adipocytes",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:21951526",
        },
        {
            "gene": "adipoq",
            "alias": "adiponectin",
            "role": "Adipocyte endocrine marker (adipokine); secreted by mature adipocytes",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:21951526",
        },
        {
            "gene": "lep",
            "alias": "leptin",
            "role": "Adipocyte endocrine factor; energy balance regulation",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:21951526",
        },
        {
            "gene": "cfd",
            "alias": "complement factor D / adipsin",
            "role": "Mature adipocyte secreted factor",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:21951526",
        },
        {
            "gene": "pnpla2",
            "alias": "ATGL (adipose triglyceride lipase)",
            "role": "Adipocyte lipid droplet hydrolysis; lipolysis marker",
            "transgenic": "ISH / qPCR",
            "refs": "PMID:28479879",
        },
    ],
}

# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def get_marker_genes_for_tissue(tissue):
    """Return a flat list of gene symbols for a given tissue."""
    entries = ZEBRAFISH_TISSUE_MARKERS.get(tissue, [])
    return [e["gene"] for e in entries]


def get_all_marker_genes():
    """Return {tissue: [gene, ...]} for every tissue."""
    return {
        tissue: [e["gene"] for e in entries]
        for tissue, entries in ZEBRAFISH_TISSUE_MARKERS.items()
    }


def get_flat_marker_set():
    """Return the union of all marker gene symbols across tissues."""
    genes = set()
    for entries in ZEBRAFISH_TISSUE_MARKERS.values():
        for e in entries:
            genes.add(e["gene"])
    return genes


def print_summary() -> None:
    """Print a quick summary table."""
    print(f"{'Tissue':<30s} {'# Markers':>10s}")
    print("-" * 42)
    total = 0
    for tissue, entries in ZEBRAFISH_TISSUE_MARKERS.items():
        print(f"{tissue:<30s} {len(entries):>10d}")
        total += len(entries)
    print("-" * 42)
    print(f"{'TOTAL':<30s} {total:>10d}")
    print(f"Unique genes: {len(get_flat_marker_set())}")


if __name__ == "__main__":
    print_summary()
    print()
    # Quick example: list all cardiovascular markers
    for tissue, genes in get_all_marker_genes().items():
        print(f"{tissue}: {', '.join(genes)}")
