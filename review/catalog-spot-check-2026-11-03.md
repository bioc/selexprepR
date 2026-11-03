# Verifica manuale a campione del catalogo selexprepR

Data di estrazione: 2026-11-03. Dieci record estratti casualmente senza
ripetizioni dai 240 record del CSV distribuito con il pacchetto.
Non sono stati esclusi record in base alla classe, alla disponibilità delle
fonti o allo stato di curation. Questo è l’unico campione estratto.

Stato: **da verificare dal maintainer**. I valori e le indicazioni delle
fonti qui sotto sono trascritti dal catalogo e dalla sua provenienza;
non sono risultati di una nuova verifica delle fonti.

Seed riproducibile: `5907265769135269765`.
Metodo: `random.Random(int(seed)).sample(rows, 10)` nell’ordine originale del CSV.
CSV SHA-256: `121b599786a718d877bd950660a3a4e0ec178db2f1b54636815d575185c3e0fe`.
Provenienza dei riferimenti: Python selexprep v0.4.6, commit
`63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`,
`src/selexprep/catalog/data/curated_metadata.json`.

Per ogni campo, annotare uno tra: **confermato**, **discrepanza**,
**non verificabile** e aggiungere la sezione/pagina o una breve nota.
Un valore vuoto con `not_stated` indica informazione non riportata nelle
fonti usate per la curation: non equivale a zero o ad assenza sperimentale.
La concordanza fra agenti non è un controllo umano.

Verificatore: ____________________
Data della verifica: ____________________

## 1. figshare:131490

Probing the Limits of Aptamer Affinity with a Microfluidic SELEX Platform

[Record originale](https://figshare.com/articles/131490)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Abstract |
| `target` | platelet derived growth factor B (PDGF-BB), thrombin, and apolipoprotein E3 (ApoE) | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Methods section |
| `target_class` | protein | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Methods section |
| `chemistry` | single-stranded DNA (ssDNA) | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Methods section |
| `n_random` | 60 randomized nucleotides | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Methods section |
| `n_rounds` | three rounds | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Methods section |
| `selection_format` | M-SELEX (microfluidic SELEX) | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0027051) — Introduction/Methods section |
| `counter_selection` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 2. PRJNA916521

Antibody-Free and Fast Profiling of Transcription Factor Binding Sites Using Aptamer-Guided Chromatin Cleavage with Sequencing (AptaChC-seq)

[Record originale](https://www.ebi.ac.uk/ena/browser/view/PRJNA916521)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | not_selex | `concordant` | [Fonte](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA916521) — BioProject title |
| `target` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `target_class` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `chemistry` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_random` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_rounds` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `selection_format` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `counter_selection` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 3. figshare:771239

Development of an Efficient Targeted Cell-SELEX Procedure for DNA Aptamer Reagents

[Record originale](https://figshare.com/articles/771239)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Results – STACS Method Overview |
| `target` | murine c-kit receptor (expressed on BJAB c-kit cells) | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Abstract – Methodology |
| `target_class` | cell surface receptor | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Discussion |
| `chemistry` | ssDNA | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Materials and Methods – Cell SELEX with DNA Aptamer Library |
| `n_random` | 29 nucleotides | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Materials and Methods – Cell SELEX with DNA Aptamer Library |
| `n_rounds` | six rounds | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Results – Cell-SELEX |
| `selection_format` | STACS (Specific TArget Cell-SELEX) | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Introduction |
| `counter_selection` | present | `concordant` | [Fonte](https://doi.org/10.1371/journal.pone.0071798) — Results – STACS Method Overview |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 4. figshare:8118461

High Affinity Aptamer for the Detection of the Biogenic Amine Histamine

[Record originale](https://figshare.com/articles/8118461)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |
| `target` | histamine | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |
| `target_class` | small molecule | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |
| `chemistry` | DNA | `concordant` | [Fonte](https://ndownloader.figshare.com/files/15135875) — Supporting Information, S1 Materials |
| `n_random` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_rounds` | ten rounds | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |
| `selection_format` | classic SELEX employing magnetic beads for target immobilization | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |
| `counter_selection` | present | `concordant` | [Fonte](https://doi.org/10.1021/acs.analchem.9b00075) — Abstract |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 5. PRJEB46967

In Vitro Selection of DNA Molecules Favouring Intramolecular Circularization by DNA Ligase

[Record originale](https://www.ebi.ac.uk/ena/browser/view/PRJEB46967)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `adjudicated` | Riferimento puntuale non registrato |
| `target` | T4DL | `single_source:codex` | [Fonte](https://doi.org/10.1021/jacs.2c12666) — Introduction |
| `target_class` | DNA ligase | `single_source:codex` | [Fonte](https://doi.org/10.1021/jacs.2c12666) — Introduction |
| `chemistry` | DNA | `concordant` | [Fonte](https://www.ebi.ac.uk/ena/browser/view/PRJEB46967) — ENA study description |
| `n_random` | N40 | `single_source:codex` | [Fonte](https://doi.org/10.1021/jacs.2c12666) — Results and Discussion, In Vitro Selection |
| `n_rounds` | 14 SELEX cycles | `single_source:codex` | [Fonte](https://doi.org/10.1021/jacs.2c12666) — Results and Discussion, In Vitro Selection |
| `selection_format` | in vitro selection | `concordant` | [Fonte](https://www.ebi.ac.uk/ena/browser/view/PRJEB46967) — ENA study description |
| `counter_selection` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 6. figshare:29756624

Supplementary file 1_Design and development of an aptamer targeting C-type lectin-like molecule-1 as a biomarker for acute myeloid leukemia: a SELEX approach.zip

[Record originale](https://figshare.com/articles/29756624)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Title |
| `target` | C-type lectin-like molecule-1 (CLL-1) | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods |
| `target_class` | protein (transmembrane glycoprotein / cell surface marker) | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods/Results |
| `chemistry` | ssDNA (76-nt library) | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods |
| `n_random` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_rounds` | eight rounds | `single_source:claude` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods |
| `selection_format` | magnetic bead-based SELEX (protein-coupled magnetic beads) | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods |
| `counter_selection` | present | `concordant` | [Fonte](https://www.frontiersin.org/articles/10.3389/fbioe.2025.1601453/full) — Methods |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 7. PRJNA957781

SELEX-seq of AP1 homodimer

[Record originale](https://www.ebi.ac.uk/ena/browser/view/PRJNA957781)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | tf_ht_selex | `concordant` | [Fonte](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957781) — BioProject description |
| `target` | AP1 (APETALA1) | `concordant` | [Fonte](https://doi.org/10.1186/s13059-025-03884-0) — Methods, SELEX-seq section |
| `target_class` | transcription factor (MADS-domain protein) | `concordant` | [Fonte](https://doi.org/10.1186/s13059-025-03884-0) — Abstract |
| `chemistry` | dsDNA | `concordant` | [Fonte](https://doi.org/10.1186/s13059-025-03884-0) — Methods, SELEX-seq section |
| `n_random` | 20 random nucleotides | `single_source:claude` | [Fonte](https://doi.org/10.1186/s13059-025-03884-0) — Methods, SELEX-seq section |
| `n_rounds` | 4 rounds | `concordant` | [Fonte](https://doi.org/10.1186/s13059-025-03884-0) — Methods, SELEX-seq section |
| `selection_format` | SELEX-seq | `adjudicated` | Riferimento puntuale non registrato |
| `counter_selection` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 8. PRJEB61115

HT-SELEX data for: GHT-SELEX demonstrates unexpectedly high intrinsic sequence specificity of human transcription factors

[Record originale](https://www.ebi.ac.uk/ena/browser/view/PRJEB61115)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | tf_ht_selex | `concordant` | [Fonte](https://www.ebi.ac.uk/ena/browser/view/PRJEB61115) — ENA Project description |
| `target` | human transcription factors (Codebook putative TF set, 331 poorly characterized putative human TFs) | `concordant` | [Fonte](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601218/) — Methods, Application of GHT-SELEX to the Codebook TF set |
| `target_class` | human transcription factors (proteins) | `concordant` | [Fonte](https://www.ebi.ac.uk/ena/browser/view/PRJEB61115) — ENA Project description |
| `chemistry` | double-stranded DNA: fragmented genomic DNA (GHT-SELEX) and random 40N synthetic ligands (HT-SELEX) | `concordant` | [Fonte](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601218/) — Methods, Application of GHT-SELEX to the Codebook TF set |
| `n_random` | random 40N ligands (HT-SELEX); GHT-SELEX uses fragmented genomic DNA with no randomized region | `concordant` | [Fonte](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601218/) — Methods, Application of GHT-SELEX to the Codebook TF set |
| `n_rounds` | number of cycles varied across 21 SELEX batches (per-experiment cycle counts recorded in Table S3) | `single_source:claude` | [Fonte](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601218/) — Methods, HT- and GHT-SELEX |
| `selection_format` | GHT-SELEX (Genomic HT-SELEX) and HT-SELEX | `concordant` | [Fonte](https://www.ebi.ac.uk/ena/browser/view/PRJEB61115) — ENA Project description |
| `counter_selection` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 9. zenodo:14272647

Structure-enhanced Deep Learning Accelerates Aptamer Selection for Small Molecule Families like Steroids

[Record originale](https://zenodo.org/records/14272647)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Methods, AptaClux: summarizing key consensus aptamer candidates from the HT-NGS result |
| `target` | hydrocortisone (cortisol) | `concordant` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Results |
| `target_class` | small molecule (steroid) | `concordant` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Abstract |
| `chemistry` | ssDNA | `single_source:claude` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Introduction |
| `n_random` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_rounds` | rounds 3, 5, and 7 (NGS sequencing rounds for this AptaVAE-guided library) | `concordant` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Methods, AptaClux: summarizing key consensus aptamer candidates from the HT-NGS result |
| `selection_format` | DL-SELEX (AptaVAE-guided SELEX with iterative rounds of qPCR, annealing, selection, and asPCR) | `concordant` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Methods, Standardized and simplified selection protocol |
| `counter_selection` | present | `single_source:claude` | [Fonte](https://doi.org/10.1093/bib/bbaf680) — Results |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________

## 10. figshare:14057938

Selection of DNA Aptamers for Root Exudate l‑Serine Using Multiple Selection Strategies

[Record originale](https://figshare.com/articles/14057938)

| Campo | Valore nel catalogo | Curation | Fonte e posizione registrate |
| --- | --- | --- | --- |
| `study_type` | aptamer_selection | `concordant` | [Fonte](https://pubmed.ncbi.nlm.nih.gov/33600189/) — Abstract |
| `target` | l-serine | `concordant` | [Fonte](https://pubmed.ncbi.nlm.nih.gov/33600189/) — Abstract |
| `target_class` | small molecule | `concordant` | [Fonte](https://pubmed.ncbi.nlm.nih.gov/33600189/) — Abstract |
| `chemistry` | DNA | `concordant` | [Fonte](https://pubmed.ncbi.nlm.nih.gov/33600189/) — Abstract |
| `n_random` | — (non riportato) | `not_stated` | Riferimento puntuale non registrato |
| `n_rounds` | Original SELEX: 11 rounds; Randomized SELEX: 8 rounds; Redirected SELEX: 5 rounds | `single_source:claude` | [Fonte](https://ndownloader.figshare.com/files/26536153) — Table S3 in Supplementary Information (jf0c06796_si_001.pdf) |
| `selection_format` | SELEX (three strategies: Original SELEX, Randomized SELEX, Redirected SELEX) | `concordant` | [Fonte](https://pubmed.ncbi.nlm.nih.gov/33600189/) — Abstract |
| `counter_selection` | present | `concordant` | [Fonte](https://ndownloader.figshare.com/files/26536153) — Table S2 and footnote in Supplementary Information (jf0c06796_si_001.pdf) |

Campi confermati: ____________________
Discrepanze o campi non verificabili: ____________________
Fonte, sezione/pagina e note: ____________________
