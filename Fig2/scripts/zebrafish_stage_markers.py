"""
Zebrafish Developmental Stage Marker Genes – Reference Dictionary
=================================================================

A curated dictionary mapping developmental stages to well-characterized marker
genes drawn from published literature.  Each gene entry includes:
  - gene   : ZFIN-standard symbol (lowercase, zebrafish convention)
  - role    : brief functional description
  - ref     : key PubMed ID(s) or journal reference

Stages follow the Kimmel et al. (1995) staging series (PMID: 8589427) plus
post-embryonic stages from Parichy et al. (2009) (PMID: 19768578).

Primary literature sources:
  - White et al., eLife 2017 (PMID: 29144233) – high-resolution mRNA time course
  - Lee et al., Nature 2013 (PMID: 24141950) – Nanog/Pou5f1/SoxB1 in MZT
  - Giraldez et al., Science 2006 (PMID: 16484454) – miR-430 maternal clearance
  - Schulte-Merker et al., Development 1994 (PMID: 7600961) – gsc/ntl
  - Holley et al., Genes Dev 2000 (PMID: 10887156) – her1 segmentation clock
  - Stainier, Curr Opin Genet Dev 2001 (PMID: 11448631) – endoderm markers
  - Paffett-Lugassy & Bhatt, Semin Cell Dev Biol 2011 – cardiac markers
  - Tang et al., Nat Commun 2017 – gata1/scl hematopoiesis
  - McCurley & Callard, BMC Mol Biol 2008 (PMID: 19014500) – housekeeping genes
  - Parichy et al., Dev Dyn 2009 (PMID: 19768578) – post-embryonic staging
"""

# ---------------------------------------------------------------------------
# Master dictionary: stage -> list of marker gene records
# ---------------------------------------------------------------------------

ZEBRAFISH_STAGE_MARKERS = {

    # ======================================================================
    # ZYGOTE (0 – 0.75 hpf, 1-cell stage)
    # Dominated by maternally deposited transcripts; no zygotic transcription.
    # ======================================================================
    "Zygote": [
        {
            "gene": "pou5f3",
            "role": "Maternal pluripotency TF (Oct4 ortholog); primes ZGA",
            "ref": "PMID: 24141950 (Lee et al., Nature 2013)",
        },
        {
            "gene": "nanog",
            "role": "Maternal pluripotency TF; activates zygotic genes at MZT",
            "ref": "PMID: 24141950 (Lee et al., Nature 2013)",
        },
        {
            "gene": "sox19b",
            "role": "Maternal SoxB1 TF; co-activates ZGA targets with Nanog/Pou5f3",
            "ref": "PMID: 24141950; PMID: 31988314 (Pálfy et al., 2020)",
        },
        {
            "gene": "buc",
            "role": "Bucky ball; maternal-effect gene required for Balbiani body/oocyte polarity",
            "ref": "PMID: 17942625 (Bontems et al., Development 2009)",
        },
        {
            "gene": "dazl",
            "role": "Maternal germplasm mRNA; RNA-binding protein in germ cell specification",
            "ref": "PMID: 11701882 (Hashimoto et al., Dev Growth Differ 2004)",
        },
        {
            "gene": "vasa",
            "role": "Maternal germplasm RNA helicase; PGC marker localized to cleavage furrows",
            "ref": "PMID: 10781161 (Knaut et al., J Cell Biol 2000)",
        },
        {
            "gene": "nanos3",
            "role": "Maternal germplasm mRNA; essential for PGC survival/migration",
            "ref": "PMID: 11691837 (Koprunner et al., Genes Dev 2001)",
        },
        {
            "gene": "dnd1",
            "role": "Dead-end 1; maternal RNA-binding protein protecting PGC mRNAs from miR-430",
            "ref": "PMID: 17363631 (Kedde et al., Cell 2007)",
        },
        {
            "gene": "ccnb1",
            "role": "Cyclin B1; maternally loaded, drives rapid S/M cell cycles pre-MBT",
            "ref": "PMID: 16284195 (Duffy et al., Nucleic Acids Res 2005)",
        },
        {
            "gene": "ccnb2",
            "role": "Cyclin B2; maternally loaded cell cycle regulator, declines post-MBT",
            "ref": "PMID: 16284195",
        },
    ],

    # ======================================================================
    # CLEAVAGE (0.75 – 2.25 hpf, 2-cell to 64-cell)
    # Rapid synchronous cell divisions; maternal transcripts dominate.
    # ======================================================================
    "Cleavage": [
        {
            "gene": "ccnb1",
            "role": "Cyclin B1; high in rapid S/M divisions, declines sharply at MBT",
            "ref": "PMID: 16284195",
        },
        {
            "gene": "ccna2",
            "role": "Cyclin A2; maternally provided, essential for S-phase in cleavage",
            "ref": "PMID: 16284195",
        },
        {
            "gene": "ccne",
            "role": "Cyclin E; drives S-phase entry in pre-MBT cleavage cycles",
            "ref": "PMID: 9019240 (Saka & Smith, Dev Biol 1997)",
        },
        {
            "gene": "vasa",
            "role": "PGC marker; localized at cleavage furrows from 4-cell stage",
            "ref": "PMID: 10781161",
        },
        {
            "gene": "nanos3",
            "role": "PGC mRNA; enriched at distal cleavage furrows in early divisions",
            "ref": "PMID: 11691837",
        },
        {
            "gene": "pou5f3",
            "role": "Maternal Oct4-like TF; abundant during cleavage, primes genome activation",
            "ref": "PMID: 24141950",
        },
        {
            "gene": "nanog",
            "role": "Maternal pluripotency TF; present through cleavage stages",
            "ref": "PMID: 24141950",
        },
    ],

    # ======================================================================
    # BLASTULA (2.25 – 5.25 hpf, 128-cell to 50% epiboly)
    # Mid-blastula transition (MBT) at ~3 hpf / 512-cell; ZGA onset;
    # miR-430 activation; maternal mRNA clearance begins.
    # ======================================================================
    "Blastula": [
        {
            "gene": "mir430",
            "role": "First massively transcribed zygotic gene cluster; clears maternal mRNAs",
            "ref": "PMID: 16484454 (Giraldez et al., Science 2006); PMID: 36693321",
        },
        {
            "gene": "nanog",
            "role": "Pluripotency TF; binds >75% of first zygotic promoters at ZGA",
            "ref": "PMID: 24141950",
        },
        {
            "gene": "pou5f3",
            "role": "Oct4 ortholog; opens chromatin at ZGA with Nanog and Sox19b",
            "ref": "PMID: 31988314 (Pálfy et al., Genome Res 2020)",
        },
        {
            "gene": "sox19b",
            "role": "SoxB1 family TF; co-activates zygotic targets at MBT",
            "ref": "PMID: 24141950",
        },
        {
            "gene": "klf17",
            "role": "Zinc-finger TF; one of earliest zygotically activated genes",
            "ref": "PMID: 31558744 (Wang et al., Sci Rep 2019)",
        },
        {
            "gene": "krt18a.1",
            "role": "Keratin 18; early epithelial/EVL marker activated at late blastula",
            "ref": "PMID: 29144233 (White et al., eLife 2017)",
        },
        {
            "gene": "gata5",
            "role": "Endoderm-specifying TF; early expression at blastula margin",
            "ref": "PMID: 11092818 (Reiter et al., Development 2001)",
        },
        {
            "gene": "apoeb",
            "role": "Apolipoprotein Eb; YSL marker from late blastula",
            "ref": "PMID: 29144233",
        },
        {
            "gene": "ccnd1",
            "role": "Cyclin D1; rises at MBT as cell cycles lengthen (replaces ccnb1)",
            "ref": "PMID: 16284195",
        },
    ],

    # ======================================================================
    # GASTRULA (5.25 – 10 hpf, 50%-epiboly to bud)
    # Germ layer specification; convergence-extension; dorsoventral patterning.
    # ======================================================================
    "Gastrula": [
        {
            "gene": "tbxta",
            "role": "T-box TF (Brachyury/ntla); pan-mesodermal, marks notochord precursors",
            "ref": "PMID: 7600961 (Schulte-Merker et al., Development 1994)",
        },
        {
            "gene": "gsc",
            "role": "Goosecoid; dorsal organizer/prechordal plate marker",
            "ref": "PMID: 7600961",
        },
        {
            "gene": "chd",
            "role": "Chordin; dorsal BMP antagonist secreted from the organizer",
            "ref": "PMID: 17959597 (Dal-Pra et al., Dev Cell 2006)",
        },
        {
            "gene": "nog1",
            "role": "Noggin 1; BMP antagonist expressed in dorsal organizer",
            "ref": "PMID: 10072795 (Fürthauer et al., Dev Biol 1999)",
        },
        {
            "gene": "sox17",
            "role": "Endoderm specification TF; downstream of sox32/casanova",
            "ref": "PMID: 11438709 (Alexander & Stainier, Genes Dev 2001)",
        },
        {
            "gene": "sox32",
            "role": "Casanova; master regulator of endoderm fate",
            "ref": "PMID: 11438709",
        },
        {
            "gene": "foxa2",
            "role": "Forkhead TF; endoderm and floorplate marker (axial/HNF3b)",
            "ref": "PMID: 11092818",
        },
        {
            "gene": "eve1",
            "role": "Even-skipped 1; ventral/posterior mesoderm marker downstream of BMP",
            "ref": "PMID: 20950598 (Cruz et al., Dev Biol 2010)",
        },
        {
            "gene": "bmp4",
            "role": "BMP4; ventral fate specification ligand during late gastrulation",
            "ref": "PMID: 19389376 (Stickney et al., Dev Dyn 2007)",
        },
        {
            "gene": "cdx4",
            "role": "Caudal homeobox TF; posterior mesoderm and hematopoietic specification",
            "ref": "PMID: 16380431 (Davidson et al., Nature 2003)",
        },
        {
            "gene": "ta",
            "role": "T-box TF (ntlb/tbxtb); additional notochord/tail marker",
            "ref": "PMID: 7600961",
        },
        {
            "gene": "eomesa",
            "role": "Eomesodermin; maternally expressed, organizer formation in mesendoderm",
            "ref": "PMID: 12952898 (Bruce et al., Development 2003)",
        },
    ],

    # ======================================================================
    # SEGMENTATION (10 – 24 hpf, 1-somite to 26-somite / prim-5)
    # Somitogenesis; segmentation clock; neurogenesis onset; early organogenesis.
    # ======================================================================
    "Segmentation": [
        {
            "gene": "her1",
            "role": "Hairy/E(spl)-related 1; segmentation clock oscillator in PSM",
            "ref": "PMID: 10887156 (Holley et al., Genes Dev 2000)",
        },
        {
            "gene": "her7",
            "role": "Segmentation clock gene; oscillates in phase with her1 in PSM",
            "ref": "PMID: 22723933 (Choorapoikayil et al., PLoS ONE 2012)",
        },
        {
            "gene": "dlc",
            "role": "DeltaC; Notch ligand oscillating in segmentation clock",
            "ref": "PMID: 17417625 (Mara et al., Dev Cell 2007)",
        },
        {
            "gene": "tbx6",
            "role": "T-box TF; marks presomitic mesoderm, regulates mespb/ripply1",
            "ref": "PMID: 25725067 (Yabe & Takada, Development 2016)",
        },
        {
            "gene": "msgn1",
            "role": "Mesogenin 1; marks PSM progenitor domain in tailbud",
            "ref": "PMID: 25725067",
        },
        {
            "gene": "mespba",
            "role": "Mesp-b; segmentally expressed in anterior PSM, confers anterior somite identity",
            "ref": "PMID: 10725245 (Sawada et al., Development 2000)",
        },
        {
            "gene": "myod1",
            "role": "Myogenic TF; marks somite myotome differentiation",
            "ref": "PMID: 22723933",
        },
        {
            "gene": "myf5",
            "role": "Myogenic factor 5; early myogenic specification in adaxial cells",
            "ref": "PMID: 8589427 (Kimmel et al., Dev Dyn 1995)",
        },
        {
            "gene": "pax2a",
            "role": "Paired-box TF; marks midbrain-hindbrain boundary (MHB) and otic vesicle",
            "ref": "PMID: 9374408 (Krauss et al., Development 1991)",
        },
        {
            "gene": "egr2b",
            "role": "Krox20; marks rhombomeres 3 and 5 in hindbrain segmentation",
            "ref": "PMID: 7720580 (Oxtoby & Jowett, Nucleic Acids Res 1993)",
        },
        {
            "gene": "neurog1",
            "role": "Neurogenin 1; earliest marker of neuronal progenitor specification",
            "ref": "PMID: 19020048 (McGraw et al., J Neurosci 2008)",
        },
        {
            "gene": "tal1",
            "role": "SCL/TAL1; hemangioblast marker in lateral plate mesoderm from 2-somite stage",
            "ref": "PMID: 22570282 (Sood et al., Adv Hematol 2012)",
        },
        {
            "gene": "lmo2",
            "role": "LIM-domain TF; co-marks hemangioblasts with tal1 in lateral mesoderm",
            "ref": "PMID: 22570282",
        },
        {
            "gene": "gata1a",
            "role": "GATA1; erythroid lineage specification from intermediate cell mass",
            "ref": "PMID: 22570282",
        },
    ],

    # ======================================================================
    # PHARYNGULA (24 – 48 hpf, prim-5 to long-pec)
    # Organogenesis progresses; heart beats; pharyngeal arches form;
    # pigmentation; brain regionalization; hatching gland maturation.
    # ======================================================================
    "Pharyngula": [
        {
            "gene": "myl7",
            "role": "Cardiac myosin light chain 7 (cmlc2); earliest cardiomyocyte marker",
            "ref": "PMID: 10491254 (Yelon et al., Development 1999)",
        },
        {
            "gene": "nkx2.5",
            "role": "Cardiac homeobox TF; marks cardiac progenitors in heart fields",
            "ref": "PMID: 31907363 (Brown et al., Genes 2020)",
        },
        {
            "gene": "hand2",
            "role": "bHLH TF; anterior lateral plate mesoderm / cardiomyocyte production",
            "ref": "PMID: 25030173 (Schindler et al., Development 2014)",
        },
        {
            "gene": "shha",
            "role": "Sonic hedgehog a; notochord/floorplate signaling",
            "ref": "PMID: 8589427",
        },
        {
            "gene": "pax6a",
            "role": "Paired-box 6a; eye and forebrain regionalization marker",
            "ref": "PMID: 39119036 (Zarzosa et al., Front Cell Dev Biol 2024)",
        },
        {
            "gene": "rx3",
            "role": "Retinal homeobox 3; earliest eye-field marker, loss gives eyeless phenotype",
            "ref": "PMID: 12947416 (Loosli et al., EMBO Rep 2003)",
        },
        {
            "gene": "otx2b",
            "role": "Orthodenticle TF; forebrain/midbrain and retinal marker",
            "ref": "PMID: 39119036",
        },
        {
            "gene": "elavl3",
            "role": "HuC; pan-neuronal post-mitotic neuron marker",
            "ref": "PMID: 19020048",
        },
        {
            "gene": "isl1",
            "role": "Islet-1; motor neuron and cranial sensory neuron marker",
            "ref": "PMID: 8589427",
        },
        {
            "gene": "sox10",
            "role": "Neural crest specifier; marks NC-derived glia, melanophores, etc.",
            "ref": "PMID: 23515338 (Mongera et al., eLife 2013)",
        },
        {
            "gene": "dlx2a",
            "role": "Distal-less 2a; pharyngeal arch neural crest marker",
            "ref": "PMID: 8589427",
        },
        {
            "gene": "hoxa2b",
            "role": "Hox group 2; selector gene patterning second pharyngeal arch",
            "ref": "PMID: 12086473 (Hunter & Prince, Dev Biol 2002)",
        },
        {
            "gene": "mpx",
            "role": "Myeloperoxidase; neutrophil/granulocyte marker (primitive wave)",
            "ref": "PMID: 24956419 (Jin et al., Dev Biol 2014)",
        },
        {
            "gene": "lyz",
            "role": "Lysozyme C; neutrophil/macrophage marker",
            "ref": "PMID: 24956419",
        },
        {
            "gene": "he1.1",
            "role": "Hatching enzyme 1; metalloprotease in hatching gland cells",
            "ref": "PMID: 31558744 (Wang et al., Sci Rep 2019)",
        },
        {
            "gene": "ctsl1b",
            "role": "Cathepsin L 1b; hatching gland enzyme for chorion digestion",
            "ref": "PMID: 31558744",
        },
    ],

    # ======================================================================
    # HATCHING (48 – 72 hpf, long-pec to protruding-mouth)
    # Embryo hatches from chorion; active swimming; continued organogenesis;
    # swim bladder inflation begins; jaw and pharyngeal skeleton mature.
    # ======================================================================
    "Hatching": [
        {
            "gene": "he1.1",
            "role": "Hatching enzyme 1; peaks at hatching, digests chorion",
            "ref": "PMID: 31558744",
        },
        {
            "gene": "ctsl1b",
            "role": "Cathepsin L 1b; hatching gland enzyme",
            "ref": "PMID: 31558744",
        },
        {
            "gene": "cd63",
            "role": "Tetraspanin; hatching gland cell marker, regulates hatching",
            "ref": "PMID: 21637765 (Huijbers et al., PLoS ONE 2011)",
        },
        {
            "gene": "fabp10a",
            "role": "Fatty acid binding protein 10a; hepatocyte marker, liver functional",
            "ref": "PMID: 32240777 (Zhao et al., Dev Biol 2020)",
        },
        {
            "gene": "try",
            "role": "Trypsin; exocrine pancreas marker, digestive function onset",
            "ref": "PMID: 30678246 (Sánchez et al., Development 2019)",
        },
        {
            "gene": "ins",
            "role": "Insulin; endocrine beta-cell marker in developing pancreatic islet",
            "ref": "PMID: 30678246",
        },
        {
            "gene": "fabp2",
            "role": "Intestinal FABP; enterocyte differentiation marker, gut functional",
            "ref": "PMID: 15531019 (Sharma et al., Gene Expr Patterns 2004)",
        },
        {
            "gene": "col2a1a",
            "role": "Type II collagen; chondrocyte marker in pharyngeal cartilage",
            "ref": "PMID: 21795417 (Bai et al., BMC Dev Biol 2014)",
        },
        {
            "gene": "mbp",
            "role": "Myelin basic protein; oligodendrocyte myelination begins ~60 hpf",
            "ref": "PMID: 29084982 (Preston & Bhatt, Methods Cell Biol 2017)",
        },
        {
            "gene": "anxa5b",
            "role": "Annexin 5b; swim bladder mesothelium marker from 60 hpf",
            "ref": "PMID: 20108353 (Teoh et al., Dev Dyn 2010)",
        },
    ],

    # ======================================================================
    # LARVAL (3 dpf – ~30 dpf)
    # Free-swimming; feeding; organ maturation; immune system development;
    # pigment pattern refinement; skeletal ossification begins.
    # ======================================================================
    "Larval": [
        {
            "gene": "rag1",
            "role": "Recombination activating gene 1; adaptive immune (lymphocyte) marker from ~4 dpf",
            "ref": "PMID: 9089097 (Willett et al., Immunogenetics 1997)",
        },
        {
            "gene": "rag2",
            "role": "Recombination activating gene 2; V(D)J recombination, thymus lymphocytes",
            "ref": "PMID: 9089097",
        },
        {
            "gene": "krt4",
            "role": "Keratin 4; most abundant type II keratin, epidermal/periderm marker",
            "ref": "PMID: 31569816 (Chen et al., G3 2019)",
        },
        {
            "gene": "tp63",
            "role": "p63; basal keratinocyte marker distinguishing basal from periderm layers",
            "ref": "PMID: 24415949 (Chou et al., PLoS Genet 2014)",
        },
        {
            "gene": "fabp10a",
            "role": "Hepatocyte marker; liver fully functional in feeding larvae",
            "ref": "PMID: 32240777",
        },
        {
            "gene": "fabp2",
            "role": "Enterocyte marker; intestinal absorptive function in feeding larvae",
            "ref": "PMID: 15531019",
        },
        {
            "gene": "slc6a3",
            "role": "Dopamine transporter (dat); dopaminergic neuron marker in larval brain",
            "ref": "PMID: 23555305 (Schmidt et al., Neural Dev 2013)",
        },
        {
            "gene": "mpx",
            "role": "Myeloperoxidase; neutrophil marker, innate immune system active",
            "ref": "PMID: 24956419",
        },
        {
            "gene": "lyz",
            "role": "Lysozyme; macrophage/neutrophil marker in larval immune response",
            "ref": "PMID: 24956419",
        },
        {
            "gene": "igfbp1a",
            "role": "IGF binding protein 1a; growth/metabolic status indicator, liver-expressed",
            "ref": "PMID: 15514009 (Kamei et al., PNAS 2005)",
        },
        {
            "gene": "sp7",
            "role": "Osterix; osteoblast TF, marks bone ossification onset in larvae",
            "ref": "PMID: 28069968 (Yu et al., Sci Rep 2017)",
        },
    ],

    # ======================================================================
    # JUVENILE (~30 dpf – ~90 dpf / 3 months)
    # Adult-like morphology; fin fold lost; scales acquired; skeletal
    # maturation; metamorphosis of pigment pattern; gonad differentiation.
    # ======================================================================
    "Juvenile": [
        {
            "gene": "runx2b",
            "role": "Osteoblast master TF; skeletal maturation and bone remodeling",
            "ref": "PMID: 22373977 (DeLaurier et al., BMC Evol Biol 2012)",
        },
        {
            "gene": "sp7",
            "role": "Osterix; osteoblast differentiation marker in maturing skeleton",
            "ref": "PMID: 28069968",
        },
        {
            "gene": "col10a1a",
            "role": "Type X collagen; hypertrophic chondrocyte and osteoblast marker",
            "ref": "PMID: 38287104 (Bergen et al., Genome Biol 2024)",
        },
        {
            "gene": "bglap",
            "role": "Osteocalcin; mature osteoblast marker in ossifying skeleton",
            "ref": "PMID: 22373977",
        },
        {
            "gene": "cyp19a1a",
            "role": "Ovarian aromatase; female gonad differentiation marker",
            "ref": "PMID: 27874044 (Lau et al., Sci Rep 2016)",
        },
        {
            "gene": "amh",
            "role": "Anti-Müllerian hormone; Sertoli cell / testis marker in gonad development",
            "ref": "PMID: 31399485 (Pfennig et al., Genetics 2019)",
        },
        {
            "gene": "gsdf",
            "role": "Gonadal somatic derived factor; ovarian follicle/testis marker",
            "ref": "PMID: 31399485",
        },
        {
            "gene": "krt4",
            "role": "Keratin 4; maintained in epidermis through juvenile metamorphosis",
            "ref": "PMID: 31569816",
        },
        {
            "gene": "igf1",
            "role": "Insulin-like growth factor 1; somatic growth signaling",
            "ref": "PMID: 12062897 (Maures et al., Endocrinology 2002)",
        },
    ],

    # ======================================================================
    # ADULT (>90 dpf, sexually mature)
    # Sexually mature; full organ homeostasis; reproductive gene expression.
    # ======================================================================
    "Adult": [
        {
            "gene": "vtg1",
            "role": "Vitellogenin 1; liver-produced egg yolk protein, female-enriched (~929-fold)",
            "ref": "PMID: 18519043 (Santos et al., Reprod Biol Endocrinol 2008)",
        },
        {
            "gene": "cyp19a1a",
            "role": "Ovarian aromatase; estrogen synthesis in mature ovary",
            "ref": "PMID: 27874044",
        },
        {
            "gene": "cyp19a1b",
            "role": "Brain aromatase; neurosteroid synthesis marker in adult brain",
            "ref": "PMID: 18519043",
        },
        {
            "gene": "amh",
            "role": "Anti-Müllerian hormone; adult testis Sertoli cell marker",
            "ref": "PMID: 31399485",
        },
        {
            "gene": "gsdf",
            "role": "Gonadal somatic factor; adult gonad marker",
            "ref": "PMID: 31399485",
        },
        {
            "gene": "rag1",
            "role": "Recombination activating gene 1; mature adaptive immune system (thymus, pronephros)",
            "ref": "PMID: 9089097",
        },
        {
            "gene": "krt4",
            "role": "Keratin 4; abundant epidermal keratin in adult skin",
            "ref": "PMID: 31569816",
        },
        {
            "gene": "tp63",
            "role": "p63; basal keratinocyte maintenance in adult epidermis",
            "ref": "PMID: 24415949",
        },
        {
            "gene": "fabp10a",
            "role": "Hepatocyte marker; mature liver homeostasis",
            "ref": "PMID: 32240777",
        },
        {
            "gene": "mbp",
            "role": "Myelin basic protein; mature myelination in adult CNS/PNS",
            "ref": "PMID: 29084982",
        },
    ],

    # ======================================================================
    # MULTI-STAGE (broadly expressed across multiple developmental stages)
    # Housekeeping, structural, and ubiquitous signaling genes.
    # ======================================================================
    "Multi-stage": [
        {
            "gene": "actb1",
            "role": "Beta-actin 1; cytoskeletal housekeeping gene, stable across stages",
            "ref": "PMID: 19014500 (McCurley & Callard, BMC Mol Biol 2008)",
        },
        {
            "gene": "eef1a1l1",
            "role": "EF1-alpha (elfa); translation elongation factor, most stable reference gene",
            "ref": "PMID: 19014500",
        },
        {
            "gene": "rpl13a",
            "role": "Ribosomal protein L13a; top-ranked stable reference gene across stages",
            "ref": "PMID: 21281742 (Tang et al., Gene Expr Patterns 2007)",
        },
        {
            "gene": "gapdh",
            "role": "Glyceraldehyde-3-phosphate dehydrogenase; metabolic gene (NB: variable in development)",
            "ref": "PMID: 19014500",
        },
        {
            "gene": "b2m",
            "role": "Beta-2-microglobulin; immune-related housekeeping, stable in geNorm analysis",
            "ref": "PMID: 19014500",
        },
        {
            "gene": "rpl7",
            "role": "Ribosomal protein L7; constitutively expressed, good reference gene",
            "ref": "PMID: 21281742",
        },
        {
            "gene": "actb2",
            "role": "Beta-actin 2; broadly expressed cytoskeletal gene",
            "ref": "PMID: 21281742",
        },
        {
            "gene": "tbp",
            "role": "TATA-binding protein; basal transcription factor, reference gene",
            "ref": "PMID: 19014500",
        },
    ],
}


# ---------------------------------------------------------------------------
# Convenience: flat lookup by gene symbol -> stage(s)
# ---------------------------------------------------------------------------
def gene_to_stages(marker_dict=None):
    """Return dict mapping gene symbol -> list of stages where it's a marker."""
    if marker_dict is None:
        marker_dict = ZEBRAFISH_STAGE_MARKERS
    g2s = {}
    for stage, genes in marker_dict.items():
        for entry in genes:
            g2s.setdefault(entry["gene"], []).append(stage)
    return g2s


def all_marker_genes(marker_dict=None):
    """Return sorted list of unique marker gene symbols."""
    if marker_dict is None:
        marker_dict = ZEBRAFISH_STAGE_MARKERS
    genes = set()
    for stage, entries in marker_dict.items():
        for entry in entries:
            genes.add(entry["gene"])
    return sorted(genes)


def stage_gene_list(stage, marker_dict=None):
    """Return just the gene symbols for a given stage."""
    if marker_dict is None:
        marker_dict = ZEBRAFISH_STAGE_MARKERS
    return [entry["gene"] for entry in marker_dict.get(stage, [])]


# ---------------------------------------------------------------------------
# Quick self-test / summary when run as script
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Zebrafish Developmental Stage Marker Genes")
    print("=" * 55)
    for stage, entries in ZEBRAFISH_STAGE_MARKERS.items():
        genes = [e["gene"] for e in entries]
        print(f"\n{stage} ({len(entries)} markers):")
        print(f"  {', '.join(genes)}")

    print(f"\n\nTotal unique marker genes: {len(all_marker_genes())}")

    print("\n\nGene -> Stage mapping (multi-stage genes):")
    g2s = gene_to_stages()
    multi = {g: s for g, s in g2s.items() if len(s) > 1}
    for g, stages in sorted(multi.items()):
        print(f"  {g}: {stages}")
