# EXTRACT_HOMOLOGS2
*Identify target homologs in amino acid datasets* 

**Disclaimer**: This repository is maintained for the express purpose of retaining bioinformatic transparency with the code used for Tassia et al. 2017.

________________________________________________________________________________________________________________________________________
## PREAMBLE ABOUT SCRIPT IN ITS ORIGINAL FORM:
The original function of this code was to identify TLR-pathway homologs across 40+ transcriptomic/genomic datasets.
The original pipeline as used for *Tassia et al. 2017* can be found as `Extract_homologues_LEGACY.sh`.

Each protein dataset is blasted against the SwissProt database and sequences which best-hit to target proteins are pulled and labeled as a putative homolog (merely by primary sequence similarity). 
To support homology given primary sequence similarity, the program annotates each sequence with `SMART` and `Pfam` domains to support homology and outputs the data in several parsable formats (`InterProScan5`'s output formats). 
Following domain annotation, the user may use user-defined domain architectures to identify complete, non-target, divergent, or fragments given knowledge/literature on known domain architectures per homolog. 
________________________________________________________________________________________________________________________________________

## ABOUT EXTRACT_HOMOLOGS2:
`EXTRACT_HOMOLOGS2` is a reimplimentation of its predacessor which was use in *Tassia et al. 2017*. Unlike its predacessor, `EXTRACT_HOMOLOGS2` has been written to be more user-friendly and (ideally) be more transportable between individual HPCs.
`EXTRACT_HOMOLOGS2` is a bash script that must be made executable in a unix environment using `chmod +x`

### DEPENDENCIES:
- [`DIAMOND`](https://github.com/bbuchfink/diamond)
- `SELECT_CONTIGS.PL` - Written by James D. White; available through this distribution
- [`CD-HIT`](https://github.com/weizhongli/cdhit)
- [`INTERPROSCAN`](https://www.ebi.ac.uk/interpro/download.html)

**Note:** If your HPC requires individual modules to be loaded into the environment, the four programs above will need to be loaded before running `EXTRACT_HOMOLOGS2`.
________________________________________________________________________________________________________________________________________
### COMMAND LINE ARGUMENTS:
- **Mandatory:**
```
  -d <str> Path to diamond database file. To make a diamond databse, use the following command: diamond makedb --in [fasta file] --db [name of databse being made]
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
________________________________________________________________________________________________________________________________________
### USAGE:

**NOTE:** For ease of use, all requisite programs should be installed into `$PATH` before using EXTRACT_HOMOLOGS2. 

**REQUIRED:** For the program to run properly, the following must be considered:
- All searchable datasets should be translated prior to running the pipeline (commonly, we use `TransDecoder` - but program preference is up to the user).
- A renaming template (`-n`) and list of target homologs (`-T`) must be provided (see examples in this distribution).
- The protein datasets must have the `.fasta` extension. **IMPORTANT:** All files in the working directory with the `.fasta` suffix will be processed by `EXTRACT_HOMOLOGS2`. Therefore, insure all `.fasta` files in the working directory are peptide files intended for processing.
- If provided, `-s` argument must terminate with `/`.

**EXAMPLE COMMAND:** 
```
EXTRACT_HOMOLOGS2.sh -d Blastdb/SwissProt.dmnd -n Naming_template.txt -T Target_homologs.txt -s Datasets/ -o TEST_OUTPUT2 -t 8 -k`
```

**PREAMBLE ON INTERPROSCAN:**
First and foremost, `Interproscan` does not impact `EXTRACT_HOMOLOGS2`'s annotation for putative homologs. Thus, if interproscan output files are not present in the `Putative_homologs/` output directory, `Interproscan` can be run on the `Putative_homologs_from_all_taxa.fasta` file after `EXTRACT_HOMOLOGS2` has finished processing.

The primary cause of `Interproscan` erroring out appears to be analysis incompatability with the host HPC. Built into the `EXTRACT_HOMOLOGS2`'s implementation of `Interproscan` are the `Pfam` and `SMART` applications (see [line 293](https://github.com/mtassia/Homolog_identification/blob/master/EXTRACT_HOMOLOGS2.sh#L293)). To ensure `EXTRACT_HOMOLOGS2` will fully complete its annotation of the identified putative homologs, I recommend running the following command to confirm command compatability *before* running `EXTRACT_HOMOLOGS2`:
```
interproscan.sh -appl SMART -appl Pfam -i [input fasta file]
```
If running this command results in an application call error (see below), try activating the `Pfam` and `SMART` application packages.
```
18/01/2019 15:09:50:456 Welcome to InterProScan-5.17-56.0
Invalid input specified for -appl/--applications parameter:
Analysis 1 does not exist or is deactivated.
```
If other errors occur when running the `interproscan` command above, please refer to the [InterProScan Wiki](https://github.com/ebi-pf-team/interproscan/wiki) and/or contact your HPC admin.

________________________________________________________________________________________________________________________________________
### OUTPUT:
**Every successful run of EXTRACT_HOMOLOGS2 (when `-k` flag isn't provided) will output the following files:**

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
