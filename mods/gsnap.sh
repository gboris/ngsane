#!/bin/bash -e

# GSNAP calling script
# author: Boris Guennewig
# date: July 2014
# modified: 

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,bwa.sh: line,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

echo ">>>>> readmapping with GSNAP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, GSNAP/GMAP index as input and 
It runs GSNAP, converts the output to .bam files, adds header information.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir

options:
  -v | --snpfile <snpfile>  SNP tolerant alignment(default: )
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
       -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group prefix
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

for MODULE in $MODULE_GSNAP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_GSNAP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BWA*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--gsnap         --\n "$(gsnap 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which gsnap)" ] && echo "[ERROR] no gsnap detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

# if [[ ! -f $FASTA/*.maps ]]; then
#     echo "[ERROR] GMAP/GSNAP index not detected. Exeute: gmap_build —D <directory> -d <genome_name> [-s none] [-k <kmer size>] <file.fasta>"
#     exit 1
# fi

# get info about input file
BAMFILE=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    BAM2BW_OPTION_ISPAIRED="True"
else
    echo "[NOTE] SINGLE library"
    BAM2BW_OPTION_ISPAIRED="False"
    PAIRED="0"
fi


# if [ -n "$READONE" ] && [ "$READONE" == "$READTWO" ]; then
#     echo "[ERROR] read1 == read2 " 1>&2 && exit 1
# elif [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ] && [ "$FORCESINGLE" = 0 ]; then
#     PAIRED="1"
#     f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
#     BAM2BW_OPTION_ISPAIRED="True"
#     echo "[NOTE] Paired library detected"
# else
#     PAIRED="0"
#     BAM2BW_OPTION_ISPAIRED="False"
#     echo "[NOTE] Single-Strand (unpaired) library detected"
# fi

## is ziped ?
ZCAT="cat" # always cat
if [[ $f = *.gz ]]; then # unless its zipped
    ZCAT="zcat" 
    ZIP="--gunzip"
else
    ZIP=""
fi

## GTF provided?
if [ -n "$GTF" ]; then
    echo "[NOTE] GTF: $GTF"
    if [ ! -f $GTF ]; then
        echo "[ERROR] GTF specified but not found!"
        exit 1
    fi 
else
    echo "[NOTE] no GTF specified!"
fi

# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi

if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi


if [ $RNA_SEQ_LIBRARY_TYPE = "fr-unstranded" ]; then
    echo "[NOTE] make bigwigs; library is fr-unstranded "
    BAM2BW_OPTION_1="FALSE"
    BAM2BW_OPTION_2="FALSE"
elif [ $RNA_SEQ_LIBRARY_TYPE = "fr-firststrand" ]; then
    echo "[NOTE] make bigwigs; library is fr-firststrand "
    BAM2BW_OPTION_1="TRUE"
    BAM2BW_OPTION_2="TRUE"
elif [ $RNA_SEQ_LIBRARY_TYPE = "fr-secondstrand" ]; then
    echo "[NOTE] make bigwigs; library is fr-secondstrand "
    BAM2BW_OPTION_1="TRUE"
    BAM2BW_OPTION_2="FALSE"	    
fi


# mkdir -p $OUTDIR

PICARD_REFERENCE=$FASTA


THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
    dmget -a $(dirname $INDEX)/*
    dmget -a $GMAP_index
    dmget -a ${f/$READONE/"*"}
	dmget -a ${OUTDIR}/$SAMPLE
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "GSNAP"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
 if [ "$PAIRED" = 1 ]; then
        echo "[NOTE] PAIRED READS"
        time gsnap \
        $ZIP \
        -D $INDEX \
        -d $GMAP_index \
        --split-output="$OUTDIR" \
        --nthreads=$CPU_GSNAP --ordered --batch=5\
        $GSNAPADDPARAM -A sam $f $f2; 
    
     else
        echo "[NOTE] SINGLE READS" 
        time gsnap \
        $ZIP \
        -D $INDEX \
        -d $GMAP_index \
        --split-output="$OUTDIR" \
        --nthreads=$CPU_GSNAP --ordered --batch=5\
        $GSNAPADDPARAM -A sam $f 
    fi
    
    # mark checkpoint
     NGSANE_CHECKPOINT_CHECK "$OUTDIR".concordant_uniq 
fi
################################################################################
NGSANE_CHECKPOINT_INIT "sort, convert into bam and merge"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

if [ "$PAIRED" = 1 ]; then
       
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".concordant_circular -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".concordant_circular
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".concordant_mult -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".concordant_mult
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".concordant_transloc -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".concordant_transloc
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".concordant_uniq -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".concordant_uniq 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".halfmapping_circular -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".halfmapping_circular 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".halfmapping_mult -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".halfmapping_mult
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".halfmapping_transloc -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".halfmapping_transloc 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".halfmapping_uniq -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".halfmapping_uniq
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".nomapping -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".nomapping 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".paired_mult -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".paired_mult 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".paired_uniq_circular -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".paired_uniq_circular 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".paired_uniq_inv -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".paired_uniq_inv
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".paired_uniq_long -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".paired_uniq_long
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".paired_uniq_scr -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".paired_uniq_scr
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_circular -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_circular
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_mult -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_mult
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_transloc -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_transloc
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_uniq -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_uniq

        
    samtools merge -@ $CPU_GSNAP "$OUTDIR".merged.bam "$OUTDIR".concordant_circular.bam "$OUTDIR".concordant_mult.bam "$OUTDIR".concordant_transloc.bam "$OUTDIR".concordant_uniq.bam "$OUTDIR".halfmapping_circular.bam "$OUTDIR".halfmapping_mult.bam "$OUTDIR".halfmapping_transloc.bam "$OUTDIR".halfmapping_uniq.bam "$OUTDIR".nomapping.bam "$OUTDIR".paired_mult.bam "$OUTDIR".paired_uniq_circular.bam "$OUTDIR".paired_uniq_inv.bam "$OUTDIR".paired_uniq_long.bam "$OUTDIR".paired_uniq_scr.bam "$OUTDIR".unpaired_circular.bam "$OUTDIR".unpaired_mult.bam "$OUTDIR".unpaired_transloc.bam "$OUTDIR".unpaired_uniq.bam
  
    
     else
         
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".nomapping -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".nomapping 
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_circular -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_circular
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_mult -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_mult
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_transloc -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_transloc
    samtools view -@ $CPU_GSNAP -bS "$OUTDIR".unpaired_uniq -o -|samtools sort -k 3,3 -k 4,4n -@ $CPU_GSNAP - "$OUTDIR".unpaired_uniq

        
    samtools merge -@ $CPU_GSNAP "$OUTDIR".merged.bam "$OUTDIR".nomapping.bam "$OUTDIR".unpaired_circular.bam "$OUTDIR".unpaired_mult.bam "$OUTDIR".unpaired_transloc.bam "$OUTDIR".unpaired_uniq.bam
    
    fi

BAMFILE="$OUTDIR".merged.bam

# rm "$OUTDIR".concordant_circular "$OUTDIR".concordant_mult "$OUTDIR".concordant_transloc "$OUTDIR".concordant_uniq "$OUTDIR".halfmapping_circular "$OUTDIR".halfmapping_mult "$OUTDIR".halfmapping_transloc "$OUTDIR".halfmapping_uniq "$OUTDIR".nomapping "$OUTDIR".paired_mult "$OUTDIR".paired_uniq_circular "$OUTDIR".paired_uniq_inv "$OUTDIR".paired_uniq_long "$OUTDIR".paired_uniq_scr "$OUTDIR".unpaired_circular "$OUTDIR".unpaired_mult "$OUTDIR".unpaired_transloc "$OUTDIR".unpaired_uniq
 
# mark checkpoint
     NGSANE_CHECKPOINT_CHECK $BAMFILE 
fi
###############################################################################
NGSANE_CHECKPOINT_INIT "clean sam"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
        INPUT=$BAMFILE \
        OUTPUT="$OUTDIR".cleaned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP

  samtools sort -@ $CPU_GSNAP "$OUTDIR".cleaned.bam "$OUTDIR"$ASD
    BAMFILE="$OUTDIR"$ASD.bam
    samtools index -@ $CPU_STAR $BAMFILE

    # mark checkpoint
     NGSANE_CHECKPOINT_CHECK $BAMFILE

fi
###############################################################################
NGSANE_CHECKPOINT_INIT "flagstat"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] samtools flagstat"
    samtools flagstat $BAMFILE > $BAMFILE.stats
    READ1=$($ZCAT $f | wc -l | gawk '{print int($1/4)}' )
    FASTQREADS=$READ1
    if [ -n "$f2" ]; then 
        READ2=$($ZCAT $f2 | wc -l | gawk '{print int($1/4)}' );
        let FASTQREADS=$READ1+$READ2
    fi
    echo $FASTQREADS" fastq reads" >> $BAMFILE.stats

    # mark checkpoint
     NGSANE_CHECKPOINT_CHECK $BAMFILE.stats

fi 
################################################################################
NGSANE_CHECKPOINT_INIT "index and calculate inner distance"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] samtools index"
    samtools index $BAMFILE

    NGSANE_CHECKPOINT_CHECK $BAMFILE.bai 

fi 
################################################################################
NGSANE_CHECKPOINT_INIT "samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] samstat"
    samstat $BAMFILE 2>&1 | tee | grep -v -P "Bad x in routine betai"
  
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $BAMFILE.stats

fi

###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

rm "$OUTDIR".concordant_circular "$OUTDIR".concordant_mult "$OUTDIR".concordant_transloc "$OUTDIR".concordant_uniq "$OUTDIR".halfmapping_circular "$OUTDIR".halfmapping_mult "$OUTDIR".halfmapping_transloc "$OUTDIR".halfmapping_uniq "$OUTDIR".nomapping "$OUTDIR".paired_mult "$OUTDIR".paired_uniq_circular "$OUTDIR".paired_uniq_inv "$OUTDIR".paired_uniq_long "$OUTDIR".paired_uniq_scr "$OUTDIR".unpaired_circular "$OUTDIR".unpaired_mult "$OUTDIR".unpaired_transloc "$OUTDIR".unpaired_uniq "$OUTDIR".cleaned.bam
 
NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE_uniq.bam.dummy ] && rm $OUTDIR/$SAMPLE_uniq.bam.dummy
echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

