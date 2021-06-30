# historical-shipping-GBS

# This repository accompanies the manuscript titled: **The reconstruction of invasion histories with genomic data in light of differing levels of anthropogenic transport** which is currently under review.

## *As such, these scripts should be treated as under review too, and are subject to change.*

- The raw fastq files produced by GBS have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB46018 (https://www.ebi.ac.uk/ena/browser/view/PRJEB46018).

- These raw reads were analysed using [ipyrad](https://ipyrad.readthedocs.io/faq.html), and the output VCF files were further filtered using [vcftools](https://vcftools.github.io/man_latest.html) as per the manuscript. 

- Final .vcf files for both species were converted to different formats for the required downstream analyses using [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/).

For R scripts for production of genomic diversity indices for both species.

For R scripts relating to the DAPC analysis.

For R scripts to reproduce the Fst heatmaps.

For R scripts for scatter plots of shipping intensity vs Fst values.
