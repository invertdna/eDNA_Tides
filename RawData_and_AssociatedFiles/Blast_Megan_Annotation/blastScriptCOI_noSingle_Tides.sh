BLAST_DB='/MBON/blastdb/nt/nt'
# BLAST PARAMETERS
PERCENT_IDENTITY="80"
WORD_SIZE="25"
EVALUE="1e-25"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="100"
CULLING="5"

################################################################################
## MEGAN ##
# For more information, see the manual provided with the software
# Specify the path to the MEGAN executable file you want to use.
# Note that in recent versions an executable was not provided; in that case, you need to reference like so: '/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub'
megan_exec='/usr/local/bin/MEGAN'

# What is the lowest taxonomic rank at which MEGAN should group OTUs?
COLLAPSE_RANK1="Family"
MINIMUM_SUPPORT="1"
MINIMUM_COMPLEXITY="0"
TOP_PERCENT="3"
MINIMUM_SUPPORT_PERCENT="0"
MINIMUM_SCORE="180"
LCA_PERCENT="80"
MAX_EXPECTED="1e-25"

# Do you want to perform a secondary MEGAN analysis, collapsing at a different taxonomic level?
PERFORM_SECONDARY_MEGAN="YES"
COLLAPSE_RANK2="Genus"


	################################################################################
	# BLAST CLUSTERS
	################################################################################
	echo $(date +%H:%M) "BLASTing..."
	blast_output="/home/mbonteam/Stanford/rpk/processed/BLASTed_COI_tides_OTUs_20170930.txt" 
	blastn \
		-query "/home/mbonteam/Stanford/rpk/raw/COI_tides_OTUs_nosingletons.fasta" \
		-db "${BLAST_DB}" \
		-num_threads 4 \
		-perc_identity "${PERCENT_IDENTITY}" \
		-word_size "${WORD_SIZE}" \
		-evalue "${EVALUE}" \
		-max_target_seqs "${MAXIMUM_MATCHES}" \
		-culling_limit="${CULLING}" \
		-outfmt "5" \
		-out "${blast_output}"


# ################################################################################
# # TAXONOMIC ANNOTATION
# ################################################################################
# 
# 
# 		echo $(date +%H:%M) 'BLAST output found; proceeding to MEGAN.'
# 		# Specify paths to megan-related files
# 		BLAST_tab="${blast_output}"
# 		MEGAN_COMMAND_FILE=/home/mbonteam/Stanford/rpk/processed/megan_commands.txt
# 		MEGAN_RMA_FILE=/home/mbonteam/Stanford/rpk/processed/meganfile.rma
# 		MEGAN_SHELL_SCRIPT=/home/mbonteam/Stanford/rpk/processed/megan_script.sh
# 
# 		echo "import blastfile='${BLAST_tab}' meganFile='${MEGAN_RMA_FILE}' \
# minScore=${MINIMUM_SCORE} \
# maxExpected=${MAX_EXPECTED} \
# topPercent=${TOP_PERCENT} \
# minSupportPercent=${MINIMUM_SUPPORT_PERCENT} \
# minSupport=${MINIMUM_SUPPORT} \
# minComplexity=${MINIMUM_COMPLEXITY} \
# lcapercent=${LCA_PERCENT} \
# blastFormat=BlastTAB;" > "${MEGAN_COMMAND_FILE}"
# 		echo "update;" >> "${MEGAN_COMMAND_FILE}"
# 		echo "collapse rank='$COLLAPSE_RANK1';" >> "${MEGAN_COMMAND_FILE}"
# 		echo "update;" >> "${MEGAN_COMMAND_FILE}"
# 		echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
# 		echo "export what=CSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK1}.csv;" >> "${MEGAN_COMMAND_FILE}"
# 		if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
# 			echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
# 			echo "update;" >> "${MEGAN_COMMAND_FILE}"
# 			echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
# 			echo "export what=CSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK2}.csv;" >> "${MEGAN_COMMAND_FILE}"
# 		fi
# 		echo "quit;" >> "${MEGAN_COMMAND_FILE}"
# 
# 		echo "#!/bin/bash" > "$MEGAN_SHELL_SCRIPT"
# 		echo "cd "${megan_exec%/*}"" >> "$MEGAN_SHELL_SCRIPT"
# 		echo "./"${megan_exec##*/}" -g -E -c ${DIR}/megan_commands.txt" >> "$MEGAN_SHELL_SCRIPT"
# 
# 		# Run MEGAN
# 		sh "${MEGAN_SHELL_SCRIPT}"

