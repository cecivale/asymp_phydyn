# Asymptomatic SARS-CoV-2 Phylodynamics

This repository contains the code and scripts used for the analyses presented in the manuscript: *Mapping Asymptomatic Transmission with Viral Phylodynamics: Insights from a Nationwide SARS-CoV-2 Screening in Young Adults*. The study investigates asymptomatic SARS-CoV-2 infections and their transmission dynamics in a cohort of Swiss men aged 20–29 with phylodynamics, using Swiss Armed Forces (SAF) screening data and community surveillance sequences.

## Repository Contents

- `asymp_saf_phydyn/` – Scripts for empirical phylodynamic analysis. 
- `asymp_simulations/` – Scripts for simulation study.  
- `README.md` – This file.

## Reproducibility

To run the Snakemake workflow, you need to include the following subworkflows in the `workflow` subfolder:  
- [talking-to-lapis](https://github.com/cecivale/talking-to-lapis)  
- [snakemake-phylo-beast2](https://github.com/cecivale/snakemake-phylo-beast2)  

## Data Availability

Publicly available sequences from the Swiss SARS-CoV-2 sequencing consortium (S3C), FOPH case data, and other demographic datasets are referenced in the manuscript.

## Contact

Cecilia Valenzuela Agüí – cecilia.valenzuela@bsse.ethz.ch
