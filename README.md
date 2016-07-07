# Homolog_identification
Identify target homologs in amino acid datasets

The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
To accomplish this task, the pipeline utilizes BLAST and the SwissProt database to identify putative homology. Following this, each putative homolog is labled with a temporary identifier denoting the blast-based homolog it best hit to. Finally, it annotates each sequence with SMART and Pfam domains to support homology and outputs the data in several formats (interproscan5's output formats). It is then on the user to visually sort homologs for domain architecture to identify complete, non-target, divergent, or fragments given knowledge on known domain architectures per homolog.
