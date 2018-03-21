#!/bin/bash


# 1- Bait minIon reads matching sequences in reference fasta file
# 2- Filter the reads using blast
# 3- Map the reads
# 4- Walk the genome to the left and the right


######################
#                    #
#    User Defined    #
#                    #
######################


# baseDir=""${HOME}"/analyses/gwalk"
baseDir=""${HOME}"/analyses/gwalk"
prefix="gwalk"

fast5="fast5_folder"
reads="gwalk_combined_trimmed.fastq.gz"
target="genomic.fasta"

#where programs are installed
prog=""${HOME}"/prog"

# Estimated genome size in bp
size=50000


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"
memCdHit=$((mem*1000))


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
logs=""${baseDir}"/logs"
basecalled=""${baseDir}"/basecalled"
qc=""${baseDir}"/qc"
matched=""${baseDir}"/matched"
trimmed=""${baseDir}"/trimmed"
assemblies=""${baseDir}"/assemblies"
mapped=""${baseDir}"/mapped"

# Create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$basecalled" ] || mkdir -p "$basecalled"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$matched" ] || mkdir -p "$matched"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$mapped" ] || mkdir -p "$mapped"


################
#              #
#     Log      #
#              #
################


#user and machine info
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

# Check for dependencies and log versions
# java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# poretools
if hash poretools 2>/dev/null; then  # if installed
    poretools -v 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "poretools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


# Canu
if hash canu 2>/dev/null; then
    canu --version | tee -a "${logs}"/log.txt
else
    echo >&2 "canu was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Racon
if hash racon 2>/dev/null; then
    echo "racon v2017-04-18" | tee -a "${logs}"/log.txt
else
    echo >&2 "racon was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# minimap
if hash minimap 2>/dev/null; then
    v=$(minimap -V)
    echo "minimap v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "minimap was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Unicycler

# Blast


# Albacore
if hash read_fast5_basecaller.py 2>/dev/null; then
    read_fast5_basecaller.py -v | tee -a "${logs}"/log.txt
else
    echo >&2 "Albacore was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi



#################
#               #
#   Functions   #
#               #
#################


# using albacore
# https://github.com/metagenomics/denbi-nanopore-training/blob/master/docs/basecalling/basecalling.rst
# https://github.com/rrwick/Basecalling-comparison

function call_bases_1D()
{

    kit="$1"  # "SQK-LSK108"
    flow="$2"  # "FLO-MIN106"
    input_fast5="$3"  # top folder - will look recursively for fast5 files because "-r" option
    output_fastq="$4"

    # will bin reads in pass/fail folder
    read_fast5_basecaller.py \
        -i "$input_fast5" \
        -t "$cpu" \
        -s "$output_fastq" \
        -k "$kit" \
        -f "$flow" \
        -r \
        -o "fastq" \
        -q 0 \
        --disable_pings

    #files of most interest are found in "${output_fastq}"/workplace/pass
}


function call_bases_1D2()
{

    kit="$1"  # "SQK-LSK108"
    flow="$2"  # "FLO-MIN106"
    input_fast5="$3"  # top folder - will look recursively for fast5 files
    output_fastq="$4"  # good reads are in"${output_fastq}"/workplace/pass

    # will bin reads in pass/fail folder
    full_1dsq_basecaller.py \
        -i "$input_fast5" \
        -t "$cpu" \
        -s "$output_fastq" \
        -k "$kit" \
        -f "$flow" \
        -r \
        -o "fastq" \
        -q 0 \
        --disable_pings
}


function demultiplex ()
{
    input_fastq="$1"
    output_folder="$2"

    [ -d "$output_folder" ] || mkdir -p "$output_folder"

    porechop     \
        -i "$input_fastq" \
        -b "$output_folder" \
        --threads "$cpu" \
        | tee "${output_folder}"_demul_porechop.log
}


function run_fastqc ()
{
    input_file="$1"
    output_folder="$2"

    #create folder for results
    [ -d "$output_folder" ] || mkdir -p "$output_folder"

    fastqc \
        --threads "$cpu" \
        --noextract \
        --o "$output_folder" \
        "$input_file"
}


function merge ()
{
    fastq_folder="$1"
    output_prefix="$2"

    mkdir -p "${fastq_folder}"/originals
    mv "${fastq_folder}"/*.fastq "${fastq_folder}"/originals

    cat "${fastq_folder}"/originals/*.fastq \
        > "${fastq_folder}"/"${output_prefix}".fastq
}


function clump ()
{
    fastq_folder="$1"  # assume ".fastq" files
    output_prefix="$2"

    #merger all fastq
    merge "$fastq_folder" "$output_prefix"

    clumpify.sh "$memJava" \
        in="${fastq_folder}"/"${output_prefix}".fastq \
        out="${fastq_folder}"/"${output_prefix}".fastq.gz \
        reorder \
        ziplevel=9

    rm "${fastq_folder}"/"${output_prefix}".fastq
}


function trim ()
{
    input_fastq="$1"

    # Trim adapters
    porechop \
        --check_reads 1000 \
        -i "$input_fastq" \
        -o "${input_fastq%.fastq.gz}_trimmed.fastq.gz" \
        --threads "$cpu" \
        | tee "${input_fastq%.fastq.gz}".porechop.log

    #remove %lines from log
    cat "${input_fastq%.fastq.gz}".porechop.log \
         | grep -v '%' \
         > "${input_fastq%.fastq.gz}".porechop.log.tmp

    mv "${input_fastq%.fastq.gz}".porechop.log.tmp \
        "${input_fastq%.fastq.gz}".porechop.log
}


function remove_chimera ()
{
    input_fastq="$1"
    output_folder="$2"

    filename="$(basename "$input_fastq")"
    name="${filename%.fastq.gz}"

    vsearch \
        --gzip_decompress \
        --threads "$cpu" \
        --log "${output_folder}"/"${name}".vsearch.log \
        --minseqlength 32 \
        --maxseqlength 50000 \
        --uchime_denovo "$input_fastq" \
        --nonchimeras "${output_folder}"/"${name}"_nonchimeras.fasta \
        --chimeras "${output_folder}"/"${name}"_chimeras.fasta
}


function bait_bbduk ()
{
    input_file="$1"  #fasta or fastq, gzipped or not
    target_fasta="$2"
    match_file="$3"
    log="$4"  # "${logs}"/bbduk_extraction.log

    bbduk.sh "$memJava" \
        overwrite=true \
        in="$input_file" \
        ref="$target_fasta" \
        threads="$cpu" \
        k=31 \
        hdist=2 \
        outm="$match_file" \
        ziplevel=9 \
        2> >(tee "$log")
}


function fastq2fasta ()
{
    input_fastq="$1"  # assume gzipped file

    zcat "$input_fastq" \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > "${input_fastq%.fastq.gz}.fasta"
}


function make_blastDB ()
{
    input_fasta="$1"

    makeblastdb \
       -in "$input_fasta" \
       -dbtype nucl \
       -parse_seqids \
       -hash_index
}


function run_blastn()
{
    query_fasta="$1"  # assume file name ends with ".fasta"
    blastDB="$2"

    blastn \
        -query "$query_fasta" \
        -db "$blastDB" \
        -out "${query_fasta%.fasta}".blastn.tsv \
        -evalue "10" \
        -word_size 28 \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore sseq' \
        -num_threads "$cpu" \
        -max_target_seqs 10 \
        -max_hsps 10


    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tevalue\tbitscore\tsseq" \
        > "${query_fasta%.fasta}".blastn.tsv.tmp
    cat "${query_fasta%.fasta}".blastn.tsv >> "${query_fasta%.fasta}".blastn.tsv.tmp
    mv "${query_fasta%.fasta}".blastn.tsv.tmp "${query_fasta%.fasta}".blastn.tsv
}


function filter_blast ()
{
    input_fasta="$1"
    target_fasta="$2"
    input_fastq="$3"  # "$reads"
    pattern="$4"  # specific pattern to extract

    run_blastn \
        "$input_fasta" \
        "$target_fasta"

    #extract matching reads
    cat "${input_fasta%.fasta}".blastn.tsv \
        | sed -e '1d' \
        | grep -F "$pattern" \
        | cut -f 1 \
        | sort | uniq \
        > "${input_fasta%.fasta}".readlist

    zcat "$input_fastq" \
        | LC_ALL=C grep --no-group-separator -A 3 -F -f "${input_fasta%.fasta}".readlist \
        | pigz > "${input_fasta%.fasta}"-filtered.fastq.gz
}


function map_graphmap()
{
    ref="$1"  # assume file ending with ".fasta"
    input_reads="$2" # assume single end gzipped fastq file or fasta

    if [ $(echo "$input_fasta" | grep -F ".fastq.gz") ]; then
        file_ext=".fastq.gz"
    elif [ $(echo "$input_fasta" | grep -F "fasta") ]; then
        file_ext=".fasta"
    else
        echo "Wrong file extension. Please use \".fasta\" or \".fastq.gz\""
        exit 1
    fi

    graphmap align \
        -t "$cpu" \
        -r "$ref" \
        -d "$input_reads" \
        -o /dev/stdout | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu" -m 10G -o "${input_reads%${file_ext}}".bam -

    #index bam file 
    samtools index "${2%${file_ext}}".bam
}


function map_bowtie2 ()
{
    ref="$1"  # assume file ending with ".fasta"
    input_reads="$2" # assume single end fastq file or fasta

    if [ $(echo "$input_reads" | grep -F ".fastq.gz") ]; then
        file_ext=".fastq.gz"
    elif [ $(echo "$input_reads" | grep -F "fasta") ]; then
        file_ext=".fasta"
    else
        echo "Wrong file extension. Please use \".fasta\" or \".fastq.gz\""
        exit 1
    fi

    # index reference
    bowtie2-build \
        "$ref" \
        "${ref%.fasta}"

    bowtie2 \
        -x "${ref%.fasta}" \
        -U "$input_reads" \
        --local \
        --very-sensitive-local \
        -p "$cpu" | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu" -m 10G -o "${input_reads%$file_ext}".bam -

    #index bam file
    samtools index "${input_reads%${file_ext}}".bam
}


function convert_bam2fastq ()
{
    input_bam="$1"
    output_fastq="$2"

    #convert bam to fastq
    bedtools bamtofastq  \
        -i "$input_bam" \
        -fq "$output_fastq"

    #compress file
    pigz -f "$output_fastq"
}


function sort_fastq_decsending ()
{
    input_fastq="$1"

    input_path=$(dirname "$input_fastq")
    input_name=$(basename "$input_fastq")

    python3.5 '/home/bioinfo/scripts/sort_fastq.py' \
    "$input_path" \
    -s "$input_name"
}


function assemble_canu ()
{
    input_fastq="$1"
    output_prefix="$2"
    output_folder="$3"
    assembly_size="$4"

    canu \
        -p "$output_prefix" \
        -d "$output_folder" \
        genomeSize="$assembly_size" \
        minReadLength=500 \
        -nanopore-raw "$input_fastq"

    # Output file is "${"$output_folder"}"/"${output_prefix}".contigs.fasta
}


function assemble_spades ()
{
    input_fastq="$1"
    output_folder="$2"

    #spades
    spades.py \
        --only-assembler \
        -k 21,33,55,77,99,127 \
        --careful \
        --iontorrent \
        -t "$cpu" \
        -o "$output_folder" \
        -s "$input_fastq"

     # Output file is "${"$output_folder"}"/contigs.fasta \
}


function assemble_unicycler ()
{
    input_fastq="$1"
    output_prefix="$2"
    output_folder="$3"
    read_type="$4"

    if [[ "$read_type" == "long" ]]; then
        opt="-l"
    elif [[ "$read_type" == "short" ]]; then
        opt="-s"
    else
        echo "Please chose \"long\" or \"short\""
        exit 1
    fi

    #Unicycler
    python3 "${prog}"/Unicycler/unicycler-runner.py \
        "$opt" "$input_fastq" \
        -o "$output_folder" \
        -t "$cpu" \
        --no_correct \
        --keep 3 \
        --verbosity 2 \
        --mode bold \
        --pilon_path "${prog}"/pilon/pilon-dev.jar

    # Output file is "${"$output_folder"}"/assembly.fasta \
}


function polish ()
{
    input_reads="$1"  # fasta format
    input_assembly="$2"
    output_folder="$3"

    filename="$(basename "$input_assembly")"
    name="${filename%.fasta}"

    map_bowtie2 "$input_assembly" "$input_reads"

    #run nanopolish
    python "${prog}"/nanopolish/scripts/nanopolish_makerange.py "$input_assembly" \
        | parallel --results "${output_folder}"/nanopolish.results -j "$cpu" \
            nanopolish variants \
                --consensus "${output_folder}"/polished.{1}.fa \
                -w {1} \
                -r "${output_folder}"/reads.fasta \
                -b "${output_folder}"/reads.sorted.bam \
                -g "$input_assembly" \
                -t 1 \
                --min-candidate-frequency 0.1

    #merge individual segments
    python "${prog}"/nanopolish/scripts/nanopolish_merge.py \
        "${output_folder}"/polished.*.fa \
        > "${output_folder}"/"${name}"_nanopolished.fasta

    #Cleanup
    rm -f "${output_folder}"/polished.*.fa
}


function run_blasr()
{
    input_fasta="$1"
    ref="$2"

    blasr \
        --nproc 48 \
        "$input_fasta" \
        "$ref" \
        > "${input_fasta%.fasta}".blasr.tsv

    echo -e "qName\ttName\tqStrand\ttStrand\tscore\tpercentSimilarity\ttStart\ttEnd\ttLength\tqStart\tqEnd\tqLength\tnCells" \
        > "${input_fasta%.fasta}".blasr.tsv.tmp
    cat "${input_fasta%.fasta}".blasr.tsv | tr " " "\t" >> "${input_fasta%.fasta}".blasr.tsv.tmp
    mv "${input_fasta%.fasta}".blasr.tsv.tmp "${input_fasta%.fasta}".blasr.tsv
}


function gwalk ()
{
    ref="$1"  # fasta files
    input_fastq="$2"  # gzipped fastq file: ".fastq.gz"
    direction="$3"
    mode="$4"

    output_folder="$(dirname "$input_fastq")"
    filename="$(basename "$input_fastq")"
    name="${filename%.fastq.gz}"

    # check that "$ref" has only one enty
    if [[ $(echo "$ref" | grep -cF ">") -gt 1 ]]; then
        echo "Please use a reference fasta file with a single entry only."
        exit 1
    fi

    # find fasta entry length
    # ref_len=$(cat "$ref" \
    #     | sed -e '1d' \
    #     | tr -d "\n" \
    #     | awk '{ print length }')

    map_bowtie2 "$ref" "$input_fastq"  # align reads

    # extract soft clipped reads on the left side of the target based on cigar
    # not in reverse complement (which are in the right side of the target)

    # $6 -> Cigar string (S=soft clipped)
    # $4 -> 1-based leftmost mapping POSition (relative to the read)
    # $2 == 16 -> reverse complement (POS 1 is still the first base on the left)
    if [ "$direction" = "left" ]; then
        bamtools convert -format sam -in "${output_folder}"/"${name}".bam \
        | awk -F $'\t' 'BEGIN {OFS = FS} {
        if($6 ~ /^[0-9][0-9][0-9]S/) {
            print
        }
        else if($6 ~ /^[0-9][0-9][0-9][0-9]S/) {
            print
        }
        else if($6 ~ /^[0-9][0-9][0-9][0-9][0-9]S/) {
            print
        }
        else if ($0 ~ /^@/) {
            print
        }
        }' | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${output_folder}"/"${name}"_"${direction}"_overhang.bam -
        samtools index "${output_folder}"/"${name}"_"${direction}"_overhang.bam

    elif [[ "$direction" == "right" ]]; then
        bamtools convert -format sam -in "${output_folder}"/"${name}".bam \
        | awk -F $'\t' 'BEGIN {OFS = FS} {
        if($6 ~ /[0-9][0-9][0-9]S$/) {
            print
        }
        else if($6 ~ /[0-9][0-9][0-9][0-9]S$/) {
            print
        }
        else if($6 ~ /[0-9][0-9][0-9][0-9][0-9]S$/) {
            print
        }
        else if ($0 ~ /^@/) {
            print
        }
        }' | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${output_folder}"/"${name}"_"${direction}"_overhang.bam -
        samtools index "${output_folder}"/"${name}"_"${direction}"_overhang.bam
    else
        echo "Please provide a direction for the walking (\"left\" or \"right\")"
        exit 1
    fi

    #visualize filtered alignment
    tablet "${output_folder}"/"${name}"_"${direction}"_overhang.bam "$target" &
}


function cluster_cd-hit-est ()
{
    input_fasta="$1"
    clustered_fasta="$2"

    #Genome wakling - reads clustering????
    # http://weizhongli-lab.org/lab-wiki/doku.php?id=cd-hit-user-guide
    # For DNAs:

        # Word size 10-11 is for thresholds 0.95 ~ 1.0
        # Word size 8,9 is for thresholds 0.90 ~ 0.95
        # Word size 7 is for thresholds 0.88 ~ 0.9
        # Word size 6 is for thresholds 0.85 ~ 0.88
        # Word size 5 is for thresholds 0.80 ~ 0.85
        # Word size 4 is for thresholds 0.75 ~ 0.8

    cd-hit-est \
        -i "$input_fasta" \
        -o "$clustered_fasta" \
        -c 0.8 \
        -s 0.8 \
        -T "$cpu" \
        -M "$memCdHit" \
        -n 4 \
        -d 0
}


############
#          #
#   Main   #
#          #
############


#starting from fast5
call_bases_1D \
    "SQK-LSK108" \
    "FLO-MIN106" \
    "$fast5" \
    "$basecalled"


# Starting from fastq
#merge and clump all fastq
clump "$read_folder" "$prefix"

#trim nanopore squencing adapter/split chimeras
trim "${read_folder}"/"${prefix}".fastq.gz

reads="${read_folder}"/"${prefix}"_trimmed.fastq.gz

# QC trimmed reads
run_fastqc \
    "$reads" \
    "$qc/fastqc/trimmed"

# Bait trimmed reads
bait_bbduk \
    "$reads" \
    "$target" \
    "${matched}"/"${prefix}".fastq.gz \
    "${matched}"/"${prefix}".log

# Make a blast database with target sequence(s)
make_blastDB "$target"

# Convert baited reads to fasta
fastq2fasta "${matched}"/"${prefix}".fastq.gz

# Filter the baited reads using blast to minimize fasle positives
# NOS, Promoter, Intron, EPSPS, 35S
filter_blast \
    "${matched}"/"${prefix}".fasta \
    "$target" \
    "$reads" \
    "NOS"

# Find reads with long left overhang
gwalk \
    "$target" \
    "${matched}"/"${prefix}"-filtered.fastq.gz \
    "left"

# Find reads with long right overhang
gwalk \
    "$target" \
    "${matched}"/"${prefix}"-filtered.fastq.gz \
    "right"

#convert bam to fastq
convert_bam2fastq \
    "${matched}"/"${prefix}"-filtered_left_overhang.bam \
    "${matched}"/"${prefix}"-filtered_left_overhang.fastq
convert_bam2fastq \
    "${matched}"/"${prefix}"-filtered_right_overhang.bam \
    "${matched}"/"${prefix}"-filtered_right_overhang.fastq

#merge both overhang fastq files
cat "${matched}"/"${prefix}"-filtered_left_overhang.fastq.gz \
    "${matched}"/"${prefix}"-filtered_right_overhang.fastq.gz \
    > "${matched}"/"${prefix}"-filtered_left_right_overhangs.fastq.gz

#assemble the reads from the genome walking
assemble_canu \
    "${matched}"/"${prefix}"-filtered_left_right_overhangs.fastq.gz \
    "$prefix" \
    "${matched}/canu" \
    "$size"

# This one seems to work the best
assemble_unicycler \
    "${matched}"/"${prefix}"-filtered_left_right_overhangs.fastq.gz \
    "$prefix" \
    "${matched}/unicycler" \
    "short"

assemble_unicycler \
    "${matched}"/"${prefix}"-filtered_left_right_overhangs.fastq.gz \
    "$prefix" \
    "${matched}/unicycler" \
    "long"

assemble_spades \
    "${matched}"/"${prefix}"-filtered_left_right_overhangs.fastq.gz \
    "${matched}/spades"
