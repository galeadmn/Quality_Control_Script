#/bin/bash
# A menu driven shell script sample template 
## ----------------------------------
# Step #1: Define variables
# ----------------------------------
RED='\033[0;41;30m'
STD='\033[0;0;39m'

FASTQC='fastqc --noextract -t 6 '
MULTIQC='multiqc '
HISAT2='hisat2 -p 85 --trim5 1 --trim3 24 '
BOWTIE2='bowtie2 -p 80 --trim5 1 --trim3 24 '

# Hisat2 indices
HISAT2_INDEX=''
HUMAN_IDX='/vol01/hisat2_indices/human/hisat2human'
RNA_IDX='/vol01/hisat2_indices/rRNA_human_mouse_rhesus/human_mouse_rhesus'
RHESUS_IDX='/vol01/hisat2_indices/rhesus/Macaca_mulatta'
MOUSE_IDX='/vol01/hisat2_indices/mouse/mouse'
RAT_IDX='/vol01/hisat2_indices/rat/hisat2rat'
ZIKA_IDX='/vol01/hisat2_indices/zika/Zika'
MASKED_IDX='/vol01/hisat2_indices/masked_rhesus/masked_rhesus'
GLOBIN_IDX='/vol01/hisat2_indices/rhesus_globin/rhesus_globins'
HCV_IDX='/vol01/hisat2_indices/HCV/HCV'
WNV_IDX='/vol01/hisat2_indices/WNV/WNV'
FLU_IDX='/vol01/hisat2_indices/FLU/FLU'
EXIT_IDX='EXIT'
INVALID_IDX='INVALID'
INDEX_NAMES=( )
INDEX_LIST=( )

# Bowtie2 indices
BT2_INDEX=''
MACAQUE_BT2_IDX='/vol01/genome/Illumina_files/macaque/Macaca_mulatta/Ensembl/Mmul_1/Sequence/Bowtie2Index/genome'
HUMAN_BT2_IDX='/vol01/genome/Illumina_files/human/Homo_sapiens/NCBI/build37.2/Sequence/Bowtie2Index/genome'
MOUSE_BT2_IDX='/vol01/genome/Illumina_files/mouse/Mus_musculus/NCBI/build37.2/Sequence/Bowtie2Index/genome'
MASKED_BT2_IDX='/vol01/genome/masked/masked'

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
	mutt -s "FastQC analysis finished for "$FASTQC_INPUT uw_galelab_hisat2output@uw.edu < /dev/null
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
	read -p "Please enter the path to the first alignment files: " ALIGN_PATH
	# Exit if there is no input entered
	if [[ -z "$ALIGN_PATH" ]] ; then
		return
	else
		while [ ! -d "$ALIGN_PATH" ] ; do
			read -p "$ALIGN_PATH does not exist. Please try again [hit return to exit]: " ALIGN_PATH
			if [[ -z $ALIGN_PATH ]] ; then
				return
			fi
		done
	fi
	format_input ALIGN_PATH
	PATH_LENGTH=${#ALIGN_PATH}

	# Prompt the user for the name of the file to record output
	read -p "Please enter the name of this dataset : " DATASET

	local counter
	counter=0
	for current_index in "${INDEX_LIST[@]}"
	do
		current_name=${INDEX_NAMES[$counter]}
		counter=$(( $counter + 1 ))

		# Create a temporary file to store intermediate results
		OUTFILE=$DATASET_$current_name.hisat2
		TEMPFILE=temp.hisat2
		echo $DATASET hisat2 results > $OUTFILE
		echo "Directory of files analyzed: " $ALIGN_PATH >> $OUTFILE
		echo "Index used for analysis : " $current_name  " " >> $OUTFILE
		echo "-------------------------------"  >> $OUTFILE
	
		# Record the start time of processing
		START=$(date +%s.%N)
		# Read the first elements of alignment paired files into an array
		firstAlignments=($ALIGN_PATH*R1_001.fastq.gz)
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

			# Remove lines reporting multiple alignment or no alignment
			sed -i '/aligned concordantly >1 times/d' $TEMPFILE
			sed -i '/aligned exactly 1 time/d' $TEMPFILE
			sed -i '/aligned >1 times/d' $TEMPFILE
			sed -i '/aligned 0 times/d' $TEMPFILE
			sed -i '/aligned concordantly 0 times/d' $TEMPFILE
			sed -i '/aligned discordantly 1 time/,+2d' $TEMPFILE
			sed -i 's/^[ \t]*//' $TEMPFILE

			# Remove the lines of just four hyphens
			sed -i '/\-\-\-\-/d' $TEMPFILE
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

process_bowtie2() {
        local choice
        local repeat_loop
        repeat_loop=1

	# Get the indices that the user wants to run bowtie2 against
        while [[ $repeat_loop -eq 1 ]]
        do
                # Prompt the user for the index file
                clear
                while [[ -z $BT2_INDEX ]]
                do
                        show_bt2_indices_menu
                        read_bt2_index_choice
                        if [[ $BT2_INDEX == $EXIT_IDX ]]
                        then
                                return
                        fi
                        if [[ $BT2_INDEX == $INVALID_IDX ]]
                        then
                                clear
                                echo "Invalid choice. Please try again."
                                BT2_INDEX=''
                        fi
                        #echo $BT2_INDEX " was chosen."
                done
                INDEX_LIST+=($BT2_INDEX)
                INDEX_NAMES+=($BT2_INDEX_NAME)
                BT2_INDEX=''
                read -p "Would you like to select another index? [y/n] " choice
                if [ ! $choice = 'y' ] && [ ! $choice = 'Y' ] ; then
                        repeat_loop=0
                fi
        done

	# Get the directory path to the fastq files to be aligned
        read -p "Please enter the path to the first alignment files: " ALIGN_PATH
        # Exit if there is no input entered
        if [[ -z "$ALIGN_PATH" ]] ; then
                return
        else
                while [ ! -d "$ALIGN_PATH" ] ; do
                        read -p "$ALIGN_PATH does not exist. Please try again [hit return to exit]: " ALIGN_PATH
                        if [[ -z $ALIGN_PATH ]] ; then
                                return
                        fi
                done
        fi
        format_input ALIGN_PATH
        PATH_LENGTH=${#ALIGN_PATH}

        # Prompt the user for the name of the file to record output
        read -p "Please enter the name of this dataset : " DATASET

        local counter
        counter=0
        for current_index in "${INDEX_LIST[@]}"
        do
                current_name=${INDEX_NAMES[$counter]}
                counter=$(( $counter + 1 ))
                # Create a temporary file to store intermediate results
                OUTFILE=$DATASET_$current_name.bowtie2
                TEMPFILE=temp.bowtie2
                echo $DATASET bowtie2 results > $OUTFILE
                echo "Index used for analysis : " $current_name  >> $OUTFILE
                echo ------------------------------   >> $OUTFILE

                # Record the start time of processing
                START=$(date +%s.%N)
                # Read the first elements of alignment paired files into an array
                firstAlignments=($ALIGN_PATH*R1_001.fastq.gz)
                for firstOfPair in "${firstAlignments[@]}"
                do
                        BaseFile=${firstOfPair%_R1_001.fastq.gz}
                        secondOfPair=$BaseFile"_R2_001.fastq.gz"
                        BaseFile=${BaseFile:$PATH_LENGTH}
                        echo >> $OUTFILE
                        echo Alignment pair $BaseFile >> $OUTFILE
         #               if [ ! -f $secondOfPair ] ; then
         #                       # Run bowtie2 on single read
         #                      bowtie2_command="$HISAT2 -x $current_index -U $firstOfPair -S /dev/null  2> $TEMPFILE"
         #               else
                                # Run bowtie2 on paired reads
                                bowtie2_command="$BOWTIE2 -x $current_index -1 $firstOfPair -2 $secondOfPair -S /dev/null  &> $TEMPFILE"
         #               fi
                        eval $bowtie2_command
         #               sed -i -e '5,15d;2,3d' $TEMPFILE
         #               sed -i -E ':a ; $!N ; s/\n\s+/ / ; ta ; P ; D' $TEMPFILE
                        cat $TEMPFILE >> $OUTFILE
                done
        END=$(date +%s.%N)
        DIFFERENCE=$(echo "$END - $START" | bc)
        echo >> $OUTFILE
        echo "Time to execute was $DIFFERENCE" >> $OUTFILE
        mutt -s $DATASET" bowtie2 analysis finished" uw_galelab_hisat2output@uw.edu < $OUTFILE
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
	echo "4. Run Bowtie2"
	echo "5. Exit"
}

# function to display hisat2 indices  menu
show_indices_menu() {
	HISAT2_INDEX=''
	echo "Select an index to use for processing:"
	echo " 1. HUMAN-MOUSE-RHESUS RNA"
	echo " 2. HUMAN"
	echo " 3. RHESUS"
	echo " 4. MOUSE"
	echo " 5. RAT"
	echo " 6. ZIKA"
	echo " 7. MASKED RHESUS"
	echo " 8. RHESUS GLOBIN"
	echo " 9. HCV"
	echo "10. WNV"
	echo "11. FLU"
	echo "12. Exit"
}

show_bt2_indices_menu() {
        BT2_INDEX=''
        echo "Select an index to use for processing:"
        echo "1. HUMAN"
        echo "2. MACAQUE"
        echo "3. MOUSE"
	echo "4. MASKED"
        echo "5. Exit"
}

read_index_choice(){
	local choice
	read -p "Enter choice [1 - 12] " choice 
	case $choice in 
		1) HISAT2_INDEX=$RNA_IDX
		   INDEX_NAME="HMR_RNA"    ;;
		2) HISAT2_INDEX=$HUMAN_IDX 
		   INDEX_NAME="HUMAN"  ;;
		3) HISAT2_INDEX=$RHESUS_IDX
		   INDEX_NAME="RHESUS" ;;
		4) HISAT2_INDEX=$MOUSE_IDX 
		   INDEX_NAME="MOUSE"  ;;
		5) HISAT2_INDEX=$RAT_IDX
		   INDEX_NAME="RAT"    ;;
		6) HISAT2_INDEX=$ZIKA_IDX
		   INDEX_NAME="ZIKA"    ;;
		7) HISAT2_INDEX=$MASKED_IDX
		   INDEX_NAME="MASKED_RHESUS"    ;;
		8) HISAT2_INDEX=$GLOBIN_IDX
		   INDEX_NAME="RHESUS_GLOBIN"    ;;
		9) HISAT2_INDEX=$HCV_IDX
		   INDEX_NAME="HCV"    ;;
		10) HISAT2_INDEX=$FLU_IDX
		   INDEX_NAME="WNV"    ;;
		11) HISAT2_INDEX=$WNV_IDX
		   INDEX_NAME="FLU"    ;;
		12) HISAT2_INDEX=$EXIT_IDX ;;
		*) exit 0 ;;
	esac
}


read_bt2_index_choice(){
        local choice
        read -p "Enter choice [1 - 5] " choice
        case $choice in
                1) BT2_INDEX=$HUMAN_BT2_IDX
                   BT2_INDEX_NAME="HUMAN"    ;;
                2) BT2_INDEX=$MACAQUE_BT2_IDX
                   BT2_INDEX_NAME="MACAQUE"  ;;
                3) BT2_INDEX=$MOUSE_BT2_IDX
                   BT2_INDEX_NAME="MOUSE"  ;;
                4) BT2_INDEX=$MASKED_BT2_IDX
                   BT2_INDEX_NAME="MASKED"  ;;
		5) BT2_INDEX=$EXIT_IDX ;;
                *) exit 0 ;;
        esac
}

# read input from the keyboard and take a action
# invoke the one() when the user select 1 from the menu option.
# invoke the two() when the user select 2 from the menu option.
# Exit when user the user select 3 form the menu option.
read_options(){
	local choice
	read -p "Enter choice [ 1 - 5] " choice
	case $choice in
		1) process_fastqc ;;
		2) process_multiqc ;;
		3) process_hisat2 ;;
		4) process_bowtie2 ;;
		5) exit 0;;
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

