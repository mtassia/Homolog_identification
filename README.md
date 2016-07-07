# Homolog_identification
Identify target homologs in amino acid datasets

The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
To accomplish this task, the pipeline utilizes BLAST and the SwissProt database to identify putative homology. Following this, each putative homolog is labled with a temporary identifier denoting the blast-based homolog it best hit to. Finally, it annotates each sequence with SMART and Pfam domains to support homology and outputs the data in several formats (interproscan5's output formats). It is then on the user to visually sort homologs for domain architecture to identify complete, non-target, divergent, or fragments given knowledge on known domain architectures per homolog.

This code is not optimized to run on other platforms, but allows for an overal understanding of the pipeline used for my own homolog identication standards. I do intend to update the code for use on other platforms.

The pipeline can be found under Extract_homologs_with_interproscan.sh. An older version of the script can be found under Extract_homologs.sh which uses HMMER to annotate Pfam domains and the SMART batch submission code - this was found to be less accurate than the interproscan verison. 

The script utilizes the following programs:
BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
SELECT_CONTIGS.PL 
CD-HIT IS INCLUDED IN THE TRANSDECODER PACKAGE AT https://github.com/TransDecoder/TransDecoder/releases
INTERPROSCAN CAN BE FOUND AT https://www.ebi.ac.uk/interpro/




Dataset_names.txt	Add files via upload	4 minutes ago
Extract_homologues.sh	Add files via upload	4 minutes ago
Extract_homologues_with_interproscan.sh	Add files via upload	4 minutes ago
Homolog_list.txt	Add files via upload	4 minutes ago
README.md	Update README.md	10 seconds ago
