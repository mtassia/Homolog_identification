# EXTRACT_HOMOLOGS2
Identify target homologs in amino acid datasets

________________________________________________________________________________________________________________________________________
## PREAMBLE ABOUT SCRIPT IN ITS ORIGINAL FORM:
The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
The original pipeline as used for *Tassia et al. 2017* can be found as `Extract_homologues_with_interproscan.sh.`

Each protein dataset is blasted against the SwissProt database and sequences which best-hit to target proteins are pulled and labeled as a putative homolog (merely by primary sequence similarity). 
To support homology given primary sequence similarity, the program annotates each sequence with SMART and Pfam domains to support homology and outputs the data in several formats (`InterProscan5`'s output formats). 
Following domain annotation, it is then on the user to visually sort homologs by domain architecture to identify complete, non-target, divergent, or fragments given knowledge on known domain architectures per homolog. 
________________________________________________________________________________________________________________________________________

## ABOUT EXTRACT_HOMOLOGS2:
`EXTRACT_HOMOLOGS2` is a reimplimentation of its predacessor which was use in *Tassia et al. 2017*. Unlike its predacessor, `EXTRACT_HOMOLOGS2` has been written to be more user-friendly and (ideally) be more transportable between individual HPCs.
`EXTRACT_HOMOLOGS2` is a bash script that must be made executable in a unix environment using `chmod +x`

### DEPENDENCIES:
- `DIAMOND` - Available at https://github.com/bbuchfink/diamond
- `SELECT_CONTIGS.PL` - Written by James D. White; available through this distribution
- `CD-HIT` - Available at https://github.com/weizhongli/cdhit
- `INTERPROSCAN` - Available at https://www.ebi.ac.uk/interpro/download.html

### COMMAND LINE ARGUMENTS:
- **Mandatory:**
```
  -d <str> Path to diamond database file
  -n <str> Naming template (example can be found in this distribution)
  -T <str> List of target SwissProt protein homologs (example can be found in this distribution)
```
- **Optional:**
```
  -h       Print help information
  -k       Keep intermediate files (disk-space expensive)
  -o <str> Output directory name (creates the directory; default value is the new directory ./Extract_homologs_output)
  -s <str> Directory containing peptide sequence files to be searched (default value is pwd)
  -S <str> Path to select_contigs.pl
  -t <int> Number of threads to be used with DIAMOND and InterproScan (default 1)
```

### USAGE:
**NOTE:** For ease of use, all requisite programs should be installed into `$PATH` before using EXTRACT_HOMOLOGS2. 

**REQUIRED:** For the program to run properly, all searchable datasets should be translated prior to running the pipeline (commonly, we use `TransDecoder` - but program preference is up to the user). Additionally, a renaming template (`-n`) and list of target homologs (`-T`) must be provided (see examples in this distribution). Finally, protein datasets must have the `.fasta` extension.

**EXAMPLE COMMAND:** 
```
Extract_homologs2.sh -d Blastdb/SwissProt.dmnd -n Naming_template.txt -T Target_homologs.txt -s Datasets/ -o TEST_OUTPUT2 -t 8 -k`
```

### OUTPUT:
**Every successful run of EXTRACT_HOMOLOGS2 will output the following files:**

- `Putative_homologs_from_all_taxa.fasta` - Sequences of putative homologs consolidated from all searched datasets
- `Putative_homologs_from_all_taxa.interpro.gff3` - InterproScan output of Putative_homologs_from_all_taxa.fasta annotation in gff3
- `Putative_homologs_from_all_taxa.interpro.tsv` - InterproScan output of Putative_homologs_from_all_taxa.fasta annotation in tsv
- `Putative_homologs_from_all_taxa.interpro.xml` - InterproScan output of Putative_homologs_from_all_taxa.fasta annotation in xml
- `Putative_homologs_from_all_taxa.blastp_vs_swissprot.outfmt6` - Diamond blastp output of Putative_homologs_from_all_taxa.fasta against SwissProt
- `Putative_homologs_from_all_taxa.besthit_vs_swissprot.outfmt6` - Diamond blastp output organized into best-hit only of Putative_homologs_from_all_taxa.fasta against SwissProt
- `[Fasta file of homologs per taxon]` - Sequences of putative homologs organized into a single fasta file per dataset
________________________________________________________________________________________________________________________________________

For any questions, suggestions, or intention of use, please contact me at mgt0007@auburn.edu

Michael Tassia,
Auburn University
