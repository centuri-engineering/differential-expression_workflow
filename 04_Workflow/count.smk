SAMPLES = expand("{samples.sample}",samples=samples.itertuples())

# ----------------------------------------------
# HISAT2-build: Indexing a reference genome
# ----------------------------------------------
rule hisat_build:
  output:
    expand( OUTPUTDIR + "{index}.1.ht2", index=index),
    expand( OUTPUTDIR + "{index}.2.ht2", index=index),
    expand( OUTPUTDIR + "{index}.3.ht2", index=index),
    expand( OUTPUTDIR + "{index}.4.ht2", index=index),
    expand( OUTPUTDIR + "{index}.5.ht2", index=index),
    expand( OUTPUTDIR + "{index}.6.ht2", index=index),
    expand( OUTPUTDIR + "{index}.7.ht2", index=index),
    expand( OUTPUTDIR + "{index}.8.ht2", index=index)

  input:
    ref = REF + config["ref"]["genome"]

  params:
    index = OUTPUTDIR + config["ref"]["index"]

  conda: 
    CONTAINER + "hisat2.yaml"

  shell:
    """
    input={input.ref}
    output={params.index}
    if [ ${{input}} == '*.gz' ];then
      gunzip ${{input}}
      input="${{input%.*}}"
    fi
    hisat2-build ${{input}} ${{output}}
    """

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
rule hisat:
  input:
    fastqc_trimmed_output= OUTPUTDIR + "03_fastqc/fastqc_trimmed_output.txt"
    
  output:
    hisat_output=OUTPUTDIR + "05_hisat/hisat.txt"

  conda: 
    CONTAINER + "hisat2.yaml"

  params:
    sample_trimmed=expand( OUTPUTDIR + "02_trimmomatic/{samples}.{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    bam = expand( OUTPUTDIR + "05_hisat/{samples}.bam", samples=SAMPLES),
    sam = expand( OUTPUTDIR + "05_hisat/{samples}.sam", samples=SAMPLES),
    reads = config["run"]["reads"],
    index = OUTPUTDIR + config["ref"]["index"]

  shell:
    """
    index={params.index}
    sample_trimmed=({params.sample_trimmed})
    len=${{#sample_trimmed[@]}}
    reads=({params.reads})
    sam=({params.sam})
    bam=({params.bam})
    flag=0
    if [ ${{reads}} == 'paired' ];then
      for (( i=0; i<$len; i=i+2 ))
        do hisat2 -p 12 -x ${{index}} -1 ${{sample_trimmed[$i]}} -2 ${{sample_trimmed[$i+1]}} -S ${{sam[$flag]}}
          samtools sort ${{sam[$flag]}} > ${{bam[$flag]}}
          rm ${{sam[$flag]}}
          flag=$((${{flag}}+1))
      done
    elif [ ${{reads}} == 'unpaired' ];then
      for (( i=0; i<$len; i++ ))
        do hisat2 -p 12 -x ${{index}} -U ${{sample_trimmed[$i]}} -S ${{sam[$i]}}
          samtools sort ${{sam[$i]}} > ${{bam[$i]}}
          rm ${{sam[$i]}    featureCounts_output = OUTPUTDIR + "05_hisat/featurecounts.txt"
}
      done
    else
      echo "Your fastq files have to be paired or unpaired. Please check the config.yaml"
    fi
    touch {output.hisat_output}
    """

# ----------------------------------------------
# featureCounts: read count/summarization program
# ----------------------------------------------

rule featureCounts:
  input:
    hisat_output = OUTPUTDIR + "05_hisat/hisat.txt"

  output:
    featureCounts_output = OUTPUTDIR + "06_featurecounts/featurecounts.txt"

  conda: 
    CONTAINER + "featureCounts.yaml"

  params:
    annotation = REF + config["ref"]["annotation"],
    bam = expand( OUTPUTDIR + "05_hisat/{samples}.bam", samples=SAMPLES),
    countmatrices = expand( OUTPUTDIR + "06_featurecounts/{samples}_count.txt", samples=SAMPLES),
    reads = config["run"]["reads"],
    geneid = config["ref"]["geneid"]

  shell:
    """
    reads=({params.reads})
    geneid=({params.geneid})
    annotation={params.annotation}
    bam=({params.bam})
    countmatrices=({params.countmatrices})
    len=${{#bam[@]}}
    mkdir 05_Output/06_featurecounts
    if [ ${{reads}} == 'paired' ];then
      for (( i=0; i<$len; i++ ))
        do featureCounts -T 12 -p -t exon -g ${{geneid}} -a ${{annotation}} -o ${{countmatrices[$i]}} ${{bam[$i]}}
      done
    elif [ ${{reads}} == 'unpaired' ];then
      for (( i=0; i<$len; i++ ))
        do featureCounts -T 12 -t exon -g ${{geneid}} -a ${{annotation}} -o ${{countmatrices[$i]}} ${{bam[$i]}}
      done
    fi
    touch {output.featureCounts_output}
    """

# ----------------------------------------------
# cpm: Counts Per Million and filtering
# ----------------------------------------------

rule cpm_filtering:
  input:
    featureCounts_output = OUTPUTDIR + "06_featurecounts/featurecounts.txt"

  output:
    cpm_filtering_output = OUTPUTDIR + "07_cpm/cpm_filtering.txt"

  params:
    count_df = report(OUTPUTDIR + "07_cpm/count.txt", caption="../report/count.rst", category="02 Count matrices"),
    output_filter_count = report(OUTPUTDIR + "07_cpm/count_filtered.txt", caption="../report/count_filtered.rst", category="02 Count matrices"),
    cpm = report(OUTPUTDIR + "07_cpm/cpm_filtered.txt", caption="../report/cpm_filtered.rst", category="02 Count matrices"),
    thresh_cpm = config["filtering"]["thresh_cpm"],
    thresh_sample = config["filtering"]["thresh_sample"],
    rmrun_list = config["filtering"]["rmrun"],
    path = OUTPUTDIR + "06_featurecounts"

  conda: 
    CONTAINER + "cpm.yaml"

  script:
    SCRIPTDIR + "cpm_filtering.R"
