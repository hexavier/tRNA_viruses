# tRNA_viruses

This repository includes the code upon submission used to reproduce the results and plots of the manuscript "Correspondence of tissue-specific translational efficiencies with the codon usage of human-infecting viruses". Data used in the code must be first downloaded in a "data" subdirectory:
- Codon usage of human-infecting viruses from CoCoPUT: https://hive.biochemistry.gwu.edu/dna.cgi?cmd=codon_usage&id=576245. Find it also in Table S1.
- List of viruses and their associtaed information: Table S1.
- Supply-to-Demand Adaptation weights of all TCGA samples: http://www.synapse.org/tRNA_TCGA
- VOGdb annotations and functional categories: http://vogdb.org/
- TCGA gene expression and clinical data: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/. In the code, this data is downloaded using the R package RTCGAToolbox (https://www.bioconductor.org/packages/release/bioc/html/RTCGAToolbox).
- tRNA expression of cells upon viral infection: Table S6.

For questions and help, please contact: xavier.hernandez@crg.eu