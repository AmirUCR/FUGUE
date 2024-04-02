# FUGUE / Fungal Universe Genome Unification Engine

FUGUE is a Python program that allows you to download genomes, CDS, GFF, and proteomes for more than 2,000 fungal species used in its parent project, [ALLEGRO](https://github.com/ucrbioinfo/allegro). It draws data from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/), [FungiDB](https://fungidb.org/), [EnsemblFungi](https://fungi.ensembl.org/index.html), and [MycoCosm](https://mycocosm.jgi.doe.gov/mycocosm/home), and removes duplicate downloads. Additionally, FUGUE can:

1. Create CDS from downloaded GFF files and mark the intron/exon boundaries
1. Utilize [DIAMOND](https://github.com/bbuchfink/diamond) to create orthogroups for your genes of interest
1. Create the required input for ALLEGRO

# Documentation
You may find the documentation for FUGUE at its [GitHub Wiki](https://github.com/AmirUCR/Fugue/wiki/).

# Support
If you run into any issues or have any suggestions for FUGUE, please report them on our GitHub Issues tracker. It's the fastest way to get support and helps us improve FUGUE for everyone.

# About
FUGUE has been developed and is maintained by Amir Mohseni, and Stefano Lonardi at the University of California, Riverside.
