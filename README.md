# EXTRACT_HOMOLOGS2
Identify target homologs in amino acid datasets

PREAMBLE ABOUT SCRIPT IN ITS ORIGINAL FORM:
The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
The original pipeline as used for Tassia et al. 2017 can be found as Extract_homologues_with_interproscan.sh.

Each protein dataset is blasted against the SwissProt database and sequences which best-hit to target proteins are pulled and labeled as a putative homolog (merely by primary sequence similarity). 
To support homology given primary sequence similarity, the program annotates each sequence with SMART and Pfam domains to support homology and outputs the data in several formats (interproscan5's output formats). 
Following domain annotation, it is then on the user to visually sort homologs by domain architecture to identify complete, non-target, divergent, or fragments given knowledge on known domain architectures per homolog. 

ABOUT EXTRACT_HOMOLOGS2:
The script utilizes the following programs:
DIAMOND - Available at https://github.com/bbuchfink/diamond
SELECT_CONTIGS.PL - Written by James D. White; available through this distribution
CD-HIT - Available at https://github.com/weizhongli/cdhit
INTERPROSCAN - Available at https://www.ebi.ac.uk/interpro/download.html

EXTRACT_HOMOLOGS2 is a reimplimentation of its predacessor which was use in Tassia et al. 2017. Unlike its predacessor, EXTRACT_HOMOLOGS2 has been written to be more user-friendly and (ideally) be more transportable between individual HPCs.
EXTRACT_HOMOLOGS2 is a bash script that must be made executable in a unix environment using 'chmod +x'

COMMAND LINE ARGUMENTS:

Mandatory:
  -d Path to diamond database file
  -n Naming template (example can be found in this distribution)
  -T List of target SwissProt protein homologs (example can be found in this distribution)
Optional:
  -h Print help information
  -k Keep intermediate files (disk-space expensive)
  -o Output directory name (creates the directory; default value is pwd)
  -s Directory containing peptide sequence files to be searched (default value is pwd)
  -S Path to select_contigs.pl
  -t Number of threads to be used with DIAMOND and InterproScan (default 1)

For ease of use, all requisite programs should be installed into $PATH before using EXTRACT_HOMOLOGS2. 

For any questions, suggestions, or intention of use, please contact me at mgt0007@auburn.edu

Michael Tassia,
Auburn University
