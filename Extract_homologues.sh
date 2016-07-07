#!/bin/bash

################################################################################################
# THIS PIPELINE IS A COMPILATION OF ALL SCRIPTS USED TO IDENTIFY THE FOLLOWING PROTIENS FROM A #
# CDS OR TRANSCRIPTOMIC DATASET.                                                               #
#                                                                       MTASSIA 7/06/2016      #
################################################################################################

# INPUT $1 = RENAMING TEMPLATE (SEE EXAMPLE INCLUDED IN DISTRIBUTION) #
# INPUT $2 = LIST OF SWISSPROT GENE IDENTIFIERS, E.G. TLR3, IRAK4, TRAF6, ETC. (SEE EXAMPLE INCLUDED IN DISTRIBUTION)#
# EXECUTE THIS SCRIPT IN A DIRECTORY CONTAINING .PEP FILES

# THIS SCRIPT UTILIZES BLASTP, SELECT_CONTIGS.PL, SMART_BATCH.PL, AND CD-HIT: 
# BLAST+ CAN BE FOUND AT https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
# SELECT_CONTIGS.PL CAN BE FOUND AT {NEEDS TO BE INCLUDED IN THE DISTRIBUTION}
# CD-HIT IS INCLUDED IN THE TRANSDECODER PACKAGE AT 

# REQUIRES SWISSPROT IN BLASTDB FORMAT - LOCATION TBD

########################################## RENAMING FUNCTION #########################################################

function rename_header {
    #PARAMETERS: 1 = RENAMING FILE (INPUT $1 ABOVE), 2 = HOMOLOGUE NAME (LINE IN INPUT $2 ABOVE), 3 = HOMOLOGUE FASTA
        
	TEMPLATE=`echo ${1}` #RENAMING TEMPLATE 
	GENE=`echo ${2}` #HOMOLOGUE TO BE RENAMED AS
	F=`echo ${3}` #FASTA FILE CONTAINING HOMOLOGOUS SEQUENCES 
        
	NEW_ID=`grep "${ID}" ${TEMPLATE} | awk '{print $2}'` # $ID COMES FROM ID VARIABLE CREATED IN MAIN
        
	awk '/^>/{print ">NEW_ID_GENE_" ++i; next}{print}' < ${F} > ${NEW_ID}_${GENE}.fasta
	sed -i "/>/s/NEW_ID/${NEW_ID}/" ${NEW_ID}_${GENE}.fasta
	sed -i "/>/s/GENE/${GENE}/" ${NEW_ID}_${GENE}.fasta
}

######################################################################################################################

# INPUT ASSIGNMENT #

RENAMING_TEMPLATE=`echo "${1}"`
HOMOLOGUE_LIST=`echo "${2}"`

# QUALITY OF LIFE CHECKS #

if [ "$#" -ne 2 ]; then #IF-STATEMENT CHECKS TO SEE IF YOU ENTERED TWO PARAMETERS
        printf 'Illegal number of parameters - two parameters must be entered:\n1) Renaming template for taxa\n2) List of genes with SwissProt identifiers\n'
        exit
fi

while read NAME #LOOP CHECKS RENAMING TEMPLATE FOR LINES WITHOUT A PAIRED OLD-NEW NAMES
do
        COLUMN_NUMBER=`echo "${NAME}" | wc -w`
        if [ "${COLUMN_NUMBER}" -ne 2 ]; then
                echo "Renaming template contains a line without a "old-new" name pair"
                exit
        fi
done < ${RENAMING_TEMPLATE}

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'` #FOR TRACKING PROGRESS OF SCRIPT
echo "${TIME_i} - Inputs entered correctly, beginning pipeline." > SCRIPT_LOG.txt

# INITIAL CLUSTERING #

for FILE in *.pep
do
	ID=`echo "${FILE}" | awk -F "." '{print $1}'` #CREATE AN ID FOR EACH PEPTIDE FILE - CORRESPONDS TO COLUMN ONE IN INPUT 1

	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'` #FOR TRACKING PROGRESS OF SCRIPT
	echo "${TIME_i} - Processing ${ID}. Beginning clustering." >> SCRIPT_LOG.txt 

	sed -i '/>/s/ .*//' ${FILE} #CLEAN HEADERS
	~/bin/cd-hit-v4.6.1-2012-08-27/cd-hit -i ${FILE} -o ${FILE}.cd-hit100 -c 1.00 #CLUSTER SEQUENCES BY 100% IDENTITY TO REMOVE FRAGMENTS OF LARGER PROTEINS
	rm -f ${FILE}.cd-hit100.clstr
	FILE=`ls | grep "${FILE}.cd-hit100"` #CHANGES FILE VARIABLE TO THE CLUSTERED DATASET
	
	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
	echo "${TIME_i} - Querying ${ID} against SwissProt database with BLASTp." >> SCRIPT_LOG.txt

# FILTER .PEP FILE BY PRIMARY SEQUENCE HOMOLOGY (BLAST) #

	blastp -query ${FILE} -db ~/bin/Swiss_Prot_DB/swissprot_db -outfmt 6 -out ${ID}.blastp_vs_swissprot.outfmt6 -evalue 1e-4 -num_threads 20 #BLAST QUERY FILE AGAINST SWISSPROT DATABASE, EVALUE CUTOFF OF 1E-4. LONGEST STEP 
	
	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
	echo "${TIME_i} - Finished blastp on ${ID}. Proceeding with BLAST output filtering" >> SCRIPT_LOG.txt 
	
	sort -k1,1 -k11,11g ${ID}.blastp_vs_swissprot.outfmt6 | sort -u -k1,1 > ${ID}.besthit_vs_swissprot.outfmt6 #PULL BEST HITS FROM THE BLAST OUTPUT BY LOWEST E-VALUE
	awk '{print $1, $2, $11}' ${ID}.besthit_vs_swissprot.outfmt6 > ${ID}.seq_hit_evalue.blast_out #REDUCE BLAST OUTPUT TO THREE COLUMNS - QUERY ID, SUBJECT ID, AND E-VALUE
	rm -f ${ID}.blastp_vs_swissprot.outfmt6
	rm -f ${ID}.besthit_vs_swissprot.outfmt6
	
# GO THROUGH THE LIST OF GENES OF INTEREST ($2) AND PULL PUTATIVE HOMOLOGUES AND RENAME FOR EACH SEQUENCE #
	
	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
	echo "${TIME_i} - Creating fasta file of homologues of interest for ${ID}" >> SCRIPT_LOG.txt
	
	while read PROT #PARSE, LINE BY LINE, A LIST OF PROTEINS HOMOLOGUES 
	do
		if [ "${PROT}" == "" ]; then
                	echo "Empty line skipped"
                	continue
        	fi
		grep "|${PROT}" ${ID}.seq_hit_evalue.blast_out | awk '{print $1}' > ${ID}.${PROT}.tmp_before_renaming #PULL QUERY SEQ HEADERS WITH BEST HIT TO PROTEINS OF INTEREST 
		select_contigs.pl -n ${ID}.${PROT}.tmp_before_renaming ${FILE} ${ID}.${PROT}.fasta.to_be_renamed #PULL CORRESPONDING SEQUENCES FROM .PEP FILE 
		rm -f ${ID}.${PROT}.tmp_before_renaming #CLEAN-UP
		rename_header ${RENAMING_TEMPLATE} ${PROT} ${ID}.${PROT}.fasta.to_be_renamed #RENAME SEQUENCES AND FILE TO INDICATE PROTEIN HOMOLOGY; ALSO CREATES $NEW_ID VARIABLE 
		rm -f ${ID}.${PROT}.fasta.to_be_renamed #CLEAN-UP
		find -empty -delete
	done < $HOMOLOGUE_LIST
	rm -f ${ID}.seq_hit_evalue.blast_out

# CONSOLIDATE HOMOLOGUE FILES INTO ONE PER TAXON #		
	
	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
	echo "${TIME_i} - Renamed ${ID} to ${NEW_ID} and compiling sequences" >> SCRIPT_LOG.txt

	ls | grep "${NEW_ID}.*fasta$" > ${NEW_ID}_homologues.txt #ASSEMBLE A FILE CONTAINING THE FILE NAMES FOR ALL PROTEINS PULLED FROM PROTEIN DATASET IN BLOCK ABOVE 
	mkdir Putative_homologues 
	while read LINE
	do
		cat $LINE >> Putative_homologues/${NEW_ID}_putative_TLR_pathway_homologues.fasta #CONSOLIDATE ALL HOMOLOGUES INTO ONE FILE PER TAXON 
	done < ${NEW_ID}_homologues.txt
	rm -f ${NEW_ID}*.fasta #CLEAN-UP
	rm -f ${NEW_ID}_homologues.txt #CLEAN-UP
	rm -f ${FILE} 
done

# NOW WORKING WITH A CONSOLIDATED ALL-TAXA DATASET #		
# ANNOTATE HOMOLOGUES WITH SMART; SMART_BATCH.PL MUST EITHER BE IN YOUR PATH OR PWD #

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Finished processing all datasets. Consolidating and annotating with SMART." >> SCRIPT_LOG.txt

cat Putative_homologues/*.fasta >> Putative_homologues_from_all_taxa.fasta #CONCATONATE ALL FASTA FILES INTO ONE
rm -rf Putative_homologues/ #CLEAN-UP
./SMART_batch.pl --includePfam --inputFile Putative_homologues_from_all_taxa.fasta #RUN REMOTE, BATCH SMART ANNOTATION ON ALL PUTATIVE HOMOLOGUES

# BLAST PUTATIVE HOMOLOGUES AGAINST SWISSPROT PRIOR TO SEQUENCE ANNOTATION TABLE CONSTRUCTION #

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Blasting consolidated dataset against SwissProt" >> SCRIPT_LOG.txt

blastp -query Putative_homologues_from_all_taxa.fasta -db ~/bin/Swiss_Prot_DB/swissprot_db -outfmt 6 -out Putative_homologues_from_all_taxa.blastp_vs_swissprot.outfmt6 -evalue 1e-4 -num_threads 20 #BLAST PUTATIVE HOMOLOGUES AGAINST SWISSPROT FOR TABLE CONSTRUCTION
sort -k1,1 -k11,11g Putative_homologues_from_all_taxa.blastp_vs_swissprot.outfmt6 | sort -u -k1,1 > Putative_homologues_from_all_taxa.besthit_vs_swissprot.outfmt6 #REDUCE BLAST OUTPUT TO BEST HITS
awk '{print $1, $2, $11}' Putative_homologues_from_all_taxa.besthit_vs_swissprot.outfmt6 > Putative_homologues_from_all_taxa.seq_hit_evalue.blast_out #REDUCE BLAST OUTPUT TO THREE COLUMNS - QUERY ID, SUBJECT ID, AND E-VALUE
rm -f Putative_homologues_from_all_taxa.blastp_vs_swissprot.outfmt6 #CLEAN-UP

# CONSTRUCT SEQUENCE ANNOTATION TABLE (CSV) #

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Compiling CSV file of putative homologues" >> SCRIPT_LOG.txt

grep '>' Putative_homologues_from_all_taxa.fasta | sed 's/^>//' > Putative_homologues_from_all_taxa.headers #MAKE A FILE CONTAINING THE LIST OF SEQUENCE HEADERS FOR PULLING BY-SEQUENCE-ANNOTATIONS
printf "ID,Best Hit vs SwissProt,E-Value of Best Hit Domain,Domain Start,Domain End,Domain E-value,Annotation Type\n" > Putative_homologues_from_all_taxa.SMART_and_BLAST_annotation.csv #CREATE HEADER FOR CSV FILE
while read LINE
do
	IDENTIFIER=`echo "${LINE}"` #USE THE LINE FROM THE HEADER LIST AS AN IDENTIFIER TO FIND SEQUENCE ANNOTATIONS IN SMART OUTPUT
	SMART_FILE=`ls SMART_results/ | grep "${IDENTIFIER}_"` #FIND THE CORRESPONDING ANNOTATION IN THE SMART OUTPUT DIRECTORY
	NUM_SMART_DOMAINS=`grep -c 'visible' SMART_results/${SMART_FILE}` #COUNT THE NUMBER OF VISIBLE DOMAINS IN THE QUERY SEQUENCE ANNOATION
	for i in $(seq 1 ${NUM_SMART_DOMAINS}) #FOR-LOOP ITERATES THROUGH EACH VISIBLE DOMAIN, NUM_SMART_DOMAINS ACTS AS UPPER BOUNDARY 
	do
		UPPER_BOUND=`grep -n 'visible' SMART_results/${SMART_FILE} | sed -n "${i}p" | awk -F":" '{print $1}'` #FIND UPPER BOUND OF i'th ANNOTATION BLOCK
        	LOWER_BOUND=`expr ${UPPER_BOUND} - 5` #FIND LOWER BOND OF i'th ANNOTATION BLOCK
        	ANNOTATION_i=`sed -n ${LOWER_BOUND},${UPPER_BOUND}p SMART_results/${SMART_FILE}` #PULL ANNOTATION BLOCK FOR THE i'th ANNOTATION IN SMART OUTPUT
        	DOMAIN=`echo "${ANNOTATION_i}" | head -1 | awk -F"=" '{print $2}'` #PULL DOMAIN TITLE
        	DOMAIN_START=`echo "${ANNOTATION_i}" | head -2 | tail -1 | awk -F"=" '{print $2}'` #PULL DOMAIN START COORDINATE
        	DOMAIN_END=`echo "${ANNOTATION_i}" | head -3 | tail -1 | awk -F"=" '{print $2}'` #PULL DOMAIN END START COORDINATE
       		DOMAIN_EVALUE=`echo "${ANNOTATION_i}" | head -4 | tail -1 | awk -F"=" '{print $2}'` #PULL DOMAIN ANNOTATION E-VALUE
        	ANNOTATION_TYPE=`echo "${ANNOTATION_i}" | tail -2 | head -1 | awk -F"=" '{print $2}'` #PULL ANNOTATION TYPE
        	BEST_HIT_VS_SWISSPROT=`grep "${IDENTIFIER} " Putative_homologues_from_all_taxa.seq_hit_evalue.blast_out | awk '{print $2}'` #PULL HOMOLOGUE'S BEST HIT FROM BLAST OUTPUT
        	BEST_HIT_VS_SWISSPROT_EVALUE=`grep "${IDENTIFIER} " Putative_homologues_from_all_taxa.seq_hit_evalue.blast_out | awk '{print $3}'` #PULL E-VALUE FOR BEST HIT IN BLAST OUTPUT
        	echo "${IDENTIFIER},${BEST_HIT_VS_SWISSPROT},${BEST_HIT_VS_SWISSPROT_EVALUE},${DOMAIN},${DOMAIN_START},${DOMAIN_END},${DOMAIN_EVALUE},${ANNOTATION_TYPE}" >> Putative_homologues_from_all_taxa.SMART_and_BLAST_annotation.csv #SAVE ANNOTATION DATA IN A CSV FILE
    	done
done < Putative_homologues_from_all_taxa.headers

# FINAL CLEAN-UP #
rm -rf SMART_results
rm -f Putative_homologues_from_all_taxa.seq_hit_evalue.blast_out
rm -f Putative_homologues_from_all_taxa.headers

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Completed pipeline. Good Job."
