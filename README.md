# st773-europe-genomic-epidemiology
Data and analysis code for the genomic epidemiology of Pseudomonas aeruginosa ST773 in the European Region, 2009–2026.
# Genomic epidemiology of Pseudomonas aeruginosa ST773 in the European Region

This repository accompanies a Dispatch by Georgios Miliotis. The study analyzed
184 ST773 genomes from 14 countries in the European Region. Recorded collection
dates spanned 2009–2026.

- `metadata/` contains the genome accession manifest and figure inputs.
- `qc/` contains assembly quality and resistance-gene summaries.
- `results/st773_parsnp_v12_184/` contains the primary phylogeny, SNP distances,
  and exploratory molecular-clock analysis.
- `results/st773_parsnp_v12_96/` contains the filtered and unfiltered SNP-distance
  matrices for the sensitivity analysis.
- `R_scripts/` and `figures/` contain the plotting code and final figures.

Public sequence records can be retrieved using the accessions in the manifest.
