# Bulbul Morphology and Range-Size Project

This repository contains the analysis workflow for a bulbul-focused project examining how morphological variation relates to species-level outcomes such as geographic range, density and niche breadth. The work brings together specimen measurements, locality data, museum holdings, collaborator datasets, phylogenies, and range/abundance information. A more detailed introduction is available in [Project Bulbul 2026.docx](Project%20Bulbul%202026.docx).

## Code files

- [code/bulbul_morphospace.Rmd](code/bulbul_morphospace.Rmd) — builds a bulbul morphospace from AVONET and related trait data and performs PCA-based exploratory analyses.
- [code/bulbul_simplecor.Rmd](code/bulbul_simplecor.Rmd) — merges range-size, abundance, and morphological datasets and prepares trait-variance summaries for downstream modelling.
- [code/bulbul_FINALCOR.R](code/bulbul_FINALCOR.R) — fits a mixed-model workflow testing whether trait coefficient of variation predicts range size while accounting for confounding effects and phylogeny.
- [code/bulbul_STAN.R](code/bulbul_STAN.R) — sets up a Stan-based phylogenetic analysis of latent variance and range-size relationships.
- [code/bulbul_sister_pair_backbone.R](code/bulbul_sister_pair_backbone.R) — provides a backbone workflow for sister-pair, range-matched analyses.
- [code/bulbul_sim1.R](code/bulbul_sim1.R) — runs a simulation study for the specimen-level trait-variance framework.
- [code/bulbul_sim2.R](code/bulbul_sim2.R) — runs a phylogenetically informed simulation study for the same framework.
- [code/duplicates.Rmd](code/duplicates.Rmd) — checks duplicate specimen identifiers and prepares duplicate-cleaning outputs.
- [code/exploratory_plotting.Rmd](code/exploratory_plotting.Rmd) — explores specimen counts, trait distributions, and related summary statistics.
- [code/Georeferencing.Rmd](code/Georeferencing.Rmd) — merges specimen data with locality information and prepares georeferenced locality outputs.
- [code/specimen_merging.Rmd](code/specimen_merging.Rmd) — combines data from collaborator museums and field projects into the main bulbul specimen dataset.
- [code/latent_variance_rangesize.stan](code/latent_variance_rangesize.stan) — defines the Stan model used for latent-variance and range-size analysis.
- [code/latent_variance_rangesize](code/latent_variance_rangesize) — compiled or intermediate Stan model artefact associated with the latent-variance analysis.
- [code/STAN_setup/STAN_main.R](code/STAN_setup/STAN_main.R) — contains the main Stan setup workflow used for model fitting.

## Data folders

- [data/](data/) — contains the main raw and processed bulbul datasets, including AVONET-derived trait tables, integrated specimen tables, locality and range-size outputs, and abundance-related files.
- [data/BIRDBASE](data/BIRDBASE) — BirdBase reference data and exported bulbul records used for taxonomic and trait cross-referencing.
- [data/Collaborator data](data/Collaborator%20data) — sample sheets and measurement files from collaborating museums, institutions, and field projects.
- [data/Density](data/Density) — bird density estimates and bulbul density datasets used for abundance-related analyses.
- [data/Duplicates](data/Duplicates) — duplicate-specimen lists and duplicate-check outputs used to clean and validate the specimen table.
- [data/Holdings](data/Holdings) — museum holding summaries, specimen-count tables, and institution-specific priority lists.
- [data/Phylogeny](data/Phylogeny) — phylogenetic trees used in the comparative analyses.
- [data/Specimen_maps](data/Specimen_maps) — spreadsheet files used to build specimen locality and museum map outputs.

## Results folders

- [results/](results/) — contains the main generated outputs from the analyses, including figures, model objects, and summary tables.
- [results/brms](results/brms) — brms model outputs and PDFs for trait-specific CV and range-size analyses.
- [results/Collaborator data comparisons](results/Collaborator%20data%20comparisons) — comparison PDFs showing how collaborator datasets align with the main specimen dataset.
- [results/Museum_map](results/Museum_map) — map-based outputs for museum and field specimen distributions.
- [results/Phylogeny](results/Phylogeny) — phylogeny-related plots, histograms, and specimen-count comparisons.
- [results/Sensitivity](results/Sensitivity) — sensitivity analyses for duplicate handling, repeatability, and sample-size effects.
- [results/Sister_pair_analysis](results/Sister_pair_analysis) — outputs and supporting files for the sister-pair analysis workflow.
- [results/Specimen_counts](results/Specimen_counts) — specimen-count figures and PDFs summarising sampling effort.
- [results/stan](results/stan) — Stan model outputs, fitted objects, and supporting model diagnostics.
