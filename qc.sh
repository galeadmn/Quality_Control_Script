#/bin/bash
# A menu driven shell script sample template 
## ----------------------------------
# Step #1: Define variables
# ----------------------------------
RED='\033[0;41;30m'
STD='\033[0;0;39m'

FASTQC='fastqc --noextract -t 6 '
MULTIQC='multiqc '
HISAT2='hisat2 -p 15 '

# Hisat2 indices
HISAT2_INDEX=''
HUMAN_IDX='/vol01/hisat2_indices/human/hisat2human'
RNA_IDX='/vol01/hisat2_indices/rRNA_human_mouse_rhesus/human_mouse_rhesus'
RHESUS_IDX='/vol01/hisat2_indices/rhesus/Macaca_mulatta'
MOUSE_IDX='/vol01/hisat2_indices/mouse/mouse'
RAT_IDX='/vol01/hisat2_indices/rat/hisat2rat'
EXIT_IDX='EXIT'
INVALID_IDX='INVALID'
INDEX_NAMES=( )
INDEX_LIST=( )


declare -a firstAlignments

# ----------------------------------
# Step #2: User defined function
# ----------------------------------
pause(){
  read -p "Press [Enter] key to continue..." fackEnterKey
}

format_input() {
  eval "last_character=${!1: -1: 1}"
  if [[ $last_character != "/" ]] ;
	then eval "$1=${!1}/" ;
  fi
}

process_fastqc(){

	# Prompt the user for the input file directory
	read -p "Please enter the path to the fastqc files: " FASTQC_INPUT
	# Exit if there is no input entered
	if [[ -z $FASTQC_INPUT ]] ; then
		return
	fi
	# Append a trailing slash if one is abscent
	format_input FASTQC_INPUT
	while [ ! -d "$FASTQC_INPUT" ] ; do
		read -p "$FASTQC_INPUT does not exist. Please try again [hit return to exit]: " FASTQC_INPUT
		if [[ -z $FASTQC_INPUT ]] ; then
			return
		fi
	done

	# Prompt the user for the output file directory
	read -p "Please enter the path to where output should be written [hit return to quit]: " FASTQC_OUTPUT
	if [[ -z $FASTQC_OUTPUT ]] ; then
		return
	fi
	format_input FASTQC_OUTPUT
	while [ ! -d "$FASTQC_OUTPUT" ] ; do
		read -p "$FASTQC_OUTPUT does not exit. Do you wish to create it?" ans
		if [[ $ans =~ ^[y,Y] ]] ; then
			mkdir -p $FASTQC_OUTPUT
		else
			read -p "Please try again [hit return to exit]: " FASTQC_OUTPUT
			if [[ -z $FASTQC_OUTPUT ]] ; then
				return
			fi
		fi
	done

	time $FASTQC -o $FASTQC_OUTPUT $FASTQC_INPUT/*.fastq.gz
	pause
}
 
process_multiqc(){

	# Prompt the user for the input file directory
	read -p "Please enter the path to the fastqc files: " MULTIQC_INPUT
	# Exit if there is no input entered
	if [[ -z $MULTIQC_INPUT ]] ; then
		return
	fi
	# Append a trailing slash if one is abscent
	format_input MULTIQC_INPUT
	while [ ! -d "$MULTIQC_INPUT" ] ; do
		read -p "$MULTIQC_INPUT does not exist. Please try again [hit return to exit]: " MULTIQC_INPUT
		if [[ -z $MULTIQC_INPUT ]] ; then
			return
		fi
	done

	# Prompt the user for the output file directory
	read -p "Please enter the path to where output should be written [hit return to quit]: " MULTIQC_OUTPUT
	if [[ -z $MULTIQC_OUTPUT ]] ; then
		return
	fi
	format_input MULTIQC_OUTPUT
	while [ ! -d "$MULTIQC_OUTPUT" ] ; do
		read -p "$MULTIQC_OUTPUT does not exit. Do you wish to create it?" ans
		if [[ $ans =~ ^[y,Y] ]] ; then
			mkdir -p $MULTIQC_OUTPUT
		else
			read -p "Please try again [hit return to exit]: " MULTIQC_OUTPUT
			if [[ -z $MULTIQC_OUTPUT ]] ; then
				return
			fi
		fi
	done

	time $MULTIQC -o $MULTIQC_OUTPUT $MULTIQC_INPUT
        pause
}

process_hisat2(){
	local choice
	local repeat_loop
	repeat_loop=1
	while [[ $repeat_loop -eq 1 ]]
	do
		# Prompt the user for the index file 
		clear
		while [[ -z $HISAT2_INDEX ]]
		do
			show_indices_menu
			read_index_choice
			if [[ $HISAT2_INDEX == $EXIT_IDX ]]
			then
				return
			fi
			if [[ $HISAT2_INDEX == $INVALID_IDX ]]
			then
				clear
				echo "Invalid choice. Please try again."
				HISAT2_INDEX=''
			fi
			#echo $HISAT2_INDEX " was chosen."
		done
		INDEX_LIST+=($HISAT2_INDEX)
		INDEX_NAMES+=($INDEX_NAME)
		HISAT2_INDEX=''
		read -p "Would you like to select another index? [y/n] " choice
		if [ ! $choice = 'y' ] && [ ! $choice = 'Y' ] ; then
			repeat_loop=0 
		fi
	done

	# Prompt the user for the directory containing alignment files 
	read -p "Please enter the path to the first alignment files: " HISAT_ALIGN_PATH
	# Exit if there is no input entered
	if [[ -z "$HISAT_ALIGN_PATH" ]] ; then
		return
	else
		while [ ! -d "$HISAT_ALIGN_PATH" ] ; do
			read -p "$HISAT_ALIGN_PATH does not exist. Please try again [hit return to exit]: " HISAT_ALIGN_PATH
			if [[ -z $HISAT_ALIGN_PATH ]] ; then
				return
			fi
		done
	fi
	format_input HISAT_ALIGN_PATH
	PATH_LENGTH=${#HISAT_ALIGN_PATH}

	# Prompt the user for the name of the file to record output
	read -p "Please enter the name of this dataset : " DATASET

	local counter
	counter=0
	for current_index in "${INDEX_LIST[@]}"
	do
		current_name=${INDEX_NAMES[$counter]}
		counter=$(( $counter + 1 ))
	#	echo $current_name
	#	read -p "press any key to continue" k
	#	continue
		# Create a temporary file to store intermediate results
		OUTFILE=$DATASET_$current_name.hisat2
		TEMPFILE=temp.hisat2
		echo $DATASET hisat2 results > $OUTFILE
		echo "Index used for analysis : " $current_name  >> $OUTFILE
		echo ------------------------------   >> $OUTFILE
	
		# Record the start time of processing
		START=$(date +%s.%N)
		# Read the first elements of alignment paired files into an array
		firstAlignments=($HISAT_ALIGN_PATH*R1_001.fastq.gz)
		for firstOfPair in "${firstAlignments[@]}" 
		do
			BaseFile=${firstOfPair%_R1_001.fastq.gz}
			secondOfPair=$BaseFile"_R2_001.fastq.gz"
			BaseFile=${BaseFile:$PATH_LENGTH}
			echo >> $OUTFILE
			echo Alignment pair $BaseFile >> $OUTFILE 
			if [ ! -f $secondOfPair ] ; then
				# Run hisat2 on single read
       		  		hisat_command="$HISAT2 -x $current_index -U $firstOfPair -S /dev/null  2> $TEMPFILE"
			else
				# Run hisat2 on paired reads
       		  		hisat_command="$HISAT2 -x $current_index -1 $firstOfPair -2 $secondOfPair -S /dev/null  2> $TEMPFILE"
			fi
			eval $hisat_command
			sed -i -e '5,15d;2,3d' $TEMPFILE
			sed -i -E ':a ; $!N ; s/\n\s+/ / ; ta ; P ; D' $TEMPFILE
			cat $TEMPFILE >> $OUTFILE
        	done
	END=$(date +%s.%N)
	DIFFERENCE=$(echo "$END - $START" | bc)
	echo >> $OUTFILE
	echo "Time to execute was $DIFFERENCE" >> $OUTFILE
	mutt -s $DATASET" hisat2 analysis finished" uw_galelab_hisat2output@uw.edu < $OUTFILE
	rm $TEMPFILE
	done
}

# function to display main menu
show_main_menu() {
	HISAT2_INDEX=''
	clear
	echo "~~~~~~~~~~~~~~~~~~~~~"	
	echo " M A I N - M E N U"
	echo "~~~~~~~~~~~~~~~~~~~~~"
	echo "1. Run FastQC"
	echo "2. Run MultiQC"
	echo "3. Run Hisat2"
	echo "4. Exit"
}

# function to display hisat2 indices  menu
show_indices_menu() {
	HISAT2_INDEX=''
	echo "Select an index to use for processing:"
	echo "1. RNA"
	echo "2. HUMAN"
	echo "3. RHESUS"
	echo "4. MOUSE"
	echo "5. RAT"
	echo "6. Exit"
}

read_index_choice(){
	local choice
	read -p "Enter choice [1 - 6] " choice 
	case $choice in 
		1) HISAT2_INDEX=$RNA_IDX
		   INDEX_NAME="RNA"    ;;
		2) HISAT2_INDEX=$HUMAN_IDX 
		   INDEX_NAME="HUMAN"  ;;
		3) HISAT2_INDEX=$RHESUS_IDX
		   INDEX_NAME="RHESUS" ;;
		4) HISAT2_INDEX=$MOUSE_IDX 
		   INDEX_NAME="MOUSE"  ;;
		5) HISAT2_INDEX=$RAT_IDX
		   INDEX_NAME="RAT"    ;;
		6) HISAT2_INDEX=$EXIT_IDX ;;
		*) exit 0 ;;
	esac
}


# read input from the keyboard and take a action
# invoke the one() when the user select 1 from the menu option.
# invoke the two() when the user select 2 from the menu option.
# Exit when user the user select 3 form the menu option.
read_options(){
	local choice
	read -p "Enter choice [ 1 - 4] " choice
	case $choice in
		1) process_fastqc ;;
		2) process_multiqc ;;
		3) process_hisat2 ;;
		4) exit 0;;
		*) echo -e "${RED}Error...${STD}" && sleep 2
	esac
}
 
# ----------------------------------------------
# Step #3: Trap CTRL+C, CTRL+Z and quit signals
# ----------------------------------------------
trap '' SIGINT SIGQUIT SIGTSTP
 
# -----------------------------------
# Step #4: Main logic - infinite loop
# ------------------------------------
while true
do
 
	show_main_menu
	read_options
done

