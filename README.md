# Detection-frontier code and data

Analysis code, derived result tables and figure source data for the manuscript

> **Uninformative non-reproduction: a pre-computable power frontier explains
> which driver–phenotype associations small cancer cohorts can and cannot test**

submitted to *BMC Medical Research Methodology*.

## What this deposit contains

A closed-form *detection frontier* `d_min(n, p)` gives the smallest standardised
effect that a two-group comparison of total size `n`, with an alteration
prevalent in a proportion `p` of the cohort, can detect at 80% power
(two-sided alpha = 0.05). The deposit contains the derivation, a simulation
validation of it under three non-normal data-generating regimes, and its
application to a complete driver-by-phenotype scan in TCGA malignant pleural
mesothelioma.

Nothing in `results/` is hand-entered. Every number traces back to a script here
or to a public data source named in the manuscript.

## Requirements

    Python >= 3.10
    numpy, scipy            # simulation validation (B7) and figures
    matplotlib              # figures only

    pip install -r requirements.txt

## Layout

    scripts/          analysis modules, run in the order below
    exploratory/      superseded analyses, NOT used in the paper (see below)
      results/        their derived tables, kept only as a record
    results/          derived tables -- every figure and table number comes from here
    figures/          vector PDF and 300-dpi PNG output
    data/             (not shipped) put downloaded public data here

## Reproducing the results

All scripts resolve paths from their own location, so the repository can be
cloned anywhere. Run from the repository root:

    python scripts/B1_power_surface.py         # power surface -> B1, B2, B3
    python scripts/B4_case_study_tests.py      # 72 driver x phenotype tests -> B4
    python scripts/B5_applicability_frontier.py  # reads B2 -> B5 tables
    python scripts/B6_figure_export.py         # reads B4 -> Figure 3 source data
    python scripts/B7_simulation_validation.py # simulation validation -> B7
    python scripts/B8_fig4_source_data.py      # needs data/ (see below)
    python scripts/F_figures.py                # regenerates Figures 1-5
    python scripts/E1_make_tables.py           # regenerates the manuscript tables

`B7` is fully self-contained and takes a few minutes. `B8` and the applied
modules in `B4` need the public data described below. `E1` and `F_figures` read
only `results/`, so both run on a fresh clone; `E1` writes its output to
`manuscript/`, which it creates if absent.

### Reproducibility of the simulation

`B7_simulation_validation.py` seeds every random draw from a CRC32 digest of the
scenario label, not from Python's built-in `hash()`. String hashing is salted per
interpreter process, so a `hash()`-derived seed changes on every run. With the
CRC32 derivation the grid regenerates bit-identically; `results/B7_*.tsv` should
be reproduced exactly.

## Data

The scripts read public data that is not redistributed here because of size.
Download into `data/`:

| Layer | Source | Expected filename |
|---|---|---|
| mRNA expression | UCSC Xena / GDC TCGA-MESO | `MESO_expr_full.tsv.gz` |
| Somatic mutation | GDC Masked Somatic Mutation | `MESO_mutations.maf` |
| Gene-level CNV | GDC ASCAT3 | `MESO_cnv_matrix.tsv.gz` |
| Protein (RPPA) | GDC | `MESO_rppa.tsv` |
| Clinical | GDC | `MESO_clinicalMatrix.tsv` |

Accession numbers and exact filters are given in the manuscript's Methods and
Availability of data section.

One module needs a directory that is **not** part of this deposit:
`B1_power_surface.py` derives the cohort sizes in `B2` from the TCGA pan-cancer
expression matrices. Point `TCGA_FULL_DIR` at a directory holding
`<COHORT>_full.tsv.gz` or `<COHORT>_HiSeqV2.gz`:

    set TCGA_FULL_DIR=D:\path\to\TCGA_full        # Windows
    export TCGA_FULL_DIR=/path/to/TCGA_full       # POSIX

If that variable is unset and no such directory exists, `B1` still writes the
power surface but leaves `results/B2_cohort_instantiation.tsv` untouched — an
empty `B2` would break `B5`, so it is never written.

## Scope note on `exploratory/`

The files in `exploratory/` belong to an earlier line of work (a
copy-number → mRNA → protein transmission analysis) that was **falsified by its
own controls** and is not reported in the manuscript. They are included only
because the negative result is itself informative and the code documents it.
`exploratory/results/` holds their derived tables for the same reason.
No result in `results/` depends on any of them. Do not cite them as part of the
paper's findings.

## Licence

Code: MIT (see `LICENSE`). Result tables and figure source data: CC BY 4.0.
