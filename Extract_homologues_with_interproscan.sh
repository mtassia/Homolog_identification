#!/bin/bash

################################################################################################
# THIS PIPELINE IS A COMPILATION OF ALL SCRIPTS USED TO IDENTIFY THE FOLLOWING PROTIENS FROM A #
# CDS OR TRANSCRIPTOMIC DATASET.                                                               #
#                                                                       MTASSIA 7/06/2016      #
################################################################################################

# INPUT $1 = RENAMING TEMPLATE (SEE EXAMPLE INCLUDED IN DISTRIBUTION) #
# INPUT $2 = LIST OF SWISSPROT GENE IDENTIFIERS, E.G. TLR3, IRAK4, TRAF6, ETC. (SEE EXAMPLE INCLUDED IN DISTRIBUTION)#
# EXECUTE THIS SCRIPT IN A DIRECTORY CONTAINING .PEP FILES

# THIS SCRIPT UTILIZES BLASTP, SELECT_CONTIGS.PL, CD-HIT, AND INTERPROSCAN: 
# BLAST+ CAN BE FOUND AT https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
# SELECT_CONTIGS.PL CAN BE FOUND AT {NEEDS TO BE INCLUDED IN THE DISTRIBUTION}
# CD-HIT IS INCLUDED IN THE TRANSDECODER PACKAGE AT https://github.com/TransDecoder/TransDecoder/releases
# INTERPROSCAN CAN BE FOUND AT https://www.ebi.ac.uk/interpro/

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
# ANNOTATE WITH INTERPROSCAN - PFAM AND SMART DATASETS #

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Finished processing all datasets. Consolidating and annotating with Interproscan." >> SCRIPT_LOG.txt

cat Putative_homologues/*.fasta >> Putative_homologues_from_all_taxa.fasta #CONCATONATE ALL FASTA FILES INTO ONE
rm -rf Putative_homologues/ #CLEAN-UP

mkdir tmp #TMP DIRECTORY FOR INTERPRO TMP FILES
WHEREAMI=`pwd` #GET THE WORKING DIRECTORY ADDRESS FOR INTERPRO - IT REQURES THE ABSOLUTE PATH
sed 's/*//g' Putative_homologues_from_all_taxa.fasta > Putative_homologues_from_all_taxa.interpro.fasta #INTERPRO REQUIRES NO "*", THIS STEP REMOVES THEM FROM THE FASTA FILE
interproscan.sh -appl SMART-6.2 -appl PfamA-27.0 -T ${WHEREAMI}/tmp -i ${WHEREAMI}/Putative_homologues_from_all_taxa.interpro.fasta

# FINAL CLEANUP #
rm -f Putative_homologues_from_all_taxa.no_stops.fasta

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Completed pipeline. Good Job."
