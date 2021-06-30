# historical-shipping-GBS
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5046438.svg)](https://doi.org/10.5281/zenodo.5046438)

# This repository accompanies the manuscript titled: **The reconstruction of invasion histories with genomic data in light of differing levels of anthropogenic transport** which is currently under review.

## *As such, these scripts should be treated as under review too, and are subject to change.*

- The raw fastq files produced by GBS have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB46018 (https://www.ebi.ac.uk/ena/browser/view/PRJEB46018).

- These raw reads were analysed using [ipyrad](https://ipyrad.readthedocs.io/faq.html), and the output VCF files were further filtered using [vcftools](https://vcftools.github.io/man_latest.html) as per the manuscript. 

- Final .vcf files for both species were converted to different formats for the required downstream analyses using [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/).

[For R scripts for production of genomic diversity indices for both species.](https://github.com/HudsonJamie/historical-shipping-GBS/tree/main/diveRsity)

[For R scripts relating to the DAPC analysis.](https://github.com/HudsonJamie/historical-shipping-GBS/tree/main/DAPC)

[For R scripts to reproduce the Fst heatmaps.](https://github.com/HudsonJamie/historical-shipping-GBS/tree/main/hierfstat)

[R scripts for scatter plots of shipping intensity vs Fst values.](https://github.com/HudsonJamie/historical-shipping-GBS/tree/main/ship_vs_genomic)
