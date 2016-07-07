# Homolog_identification
Identify target homologs in amino acid datasets

The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
To accomplish this task, the pipeline utilizes BLAST and the SwissProt database to identify putative homology. Following this, each putative homolog is labled with a temporary identifier denoting the blast-based homolog it best hit to. Finally, it annotates each sequence with SMART and Pfam domains to support homology and outputs the data in several formats (interproscan5's output formats). It is then on the user to visually sort homologs for domain architecture to identify complete, non-target, divergent, or fragments given knowledge on known domain architectures per homolog.

The pipeline can be found under Extract_homologs_with_interproscan.sh. An older version of the script can be found under Extract_homologs.sh which uses HMMER to annotate Pfam domains and the SMART batch submission code - this was found to be less accurate than the interproscan verison. 

The script utilizes the following programs:
BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
SELECT_CONTIGS.PL (included in this distribution, but this is not my own code)
CD-HIT (https://github.com/TransDecoder/TransDecoder/releases)
INTERPROSCAN (https://www.ebi.ac.uk/interpro/)

Extract_homologs_with_interproscan.sh takes two arguments:
1. Renaming template (see Dataset_names.txt as an example of this format)
2. SwissProt target ID's (See Homolog_list.txt as an example of this format)

This code is not optimized to run on other platforms, but allows for an overal understanding of the pipeline used for my own homolog identication standards. I do intend to update the code for use on other platforms.

For any questions, suggestions, or intention of use, please contact me at mgt0007@auburn.edu

Michael Tassia
Auburn University
