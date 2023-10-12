#SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
SAMPLES = expand("{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
EXT = config["run"]["ext"]

# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------

rule fastqc:
  input:
    fastqc_input= OUTPUTDIR + "01_fastqc/fastqc_input.txt"

  output:
    fastqc_output= OUTPUTDIR + "01_fastqc/fastqc_output.txt"

  params:
    rawdata=expand( "00_RawData/{samples}.{run}.{ext}", samples=SAMPLES, run=RUN, ext=EXT)

  conda: 
    CONTAINER + "fastqc.yaml"

  shell:
    """
    fastqc --outdir 05_Output/01_fastqc/ {params.rawdata}
    touch {output.fastqc_output}
    """

# ----------------------------------------------
# Trimmomatic: trimming reads and removing adapter sequences
# ----------------------------------------------

rule trimmomatic:
  input:
    fastqc_output= OUTPUTDIR + "01_fastqc/fastqc_output.txt"

  output:
    trimmomatic_output= OUTPUTDIR + "02_trimmomatic/trimmomatic_output.txt"

  params:
    sample=expand("00_RawData/{samples}.{run}.{ext}", samples=SAMPLES, run=RUN, ext=EXT),
    sample_trimmed=expand( "05_Output/02_trimmomatic/{samples}.{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    sample_untrimmed=expand( "05_Output/02_trimmomatic/{samples}.{run}un.trimmed.fastq", samples=SAMPLES, run=RUN),
    illuminaclip=config["clean"]["illuminaclip"],

  conda: 
    CONTAINER + "trimmomatic.yaml"

  shell:
    """
    sample=({params.sample})
    sample_trimmed=({params.sample_trimmed})
    sample_untrimmed=({params.sample_untrimmed})
    len=${{#sample[@]}}
    for (( i=0; i<$len; i=i+2 ))
    do trimmomatic PE -threads 4 ${{sample[$i]}} ${{sample[$i+1]}} ${{sample_trimmed[$i]}} ${{sample_untrimmed[$i]}} ${{sample_trimmed[$i+1]}} ${{sample_untrimmed[$i+1]}} ILLUMINACLIP:{params.illuminaclip}:2:30:10 LEADING:20 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
    done
    touch {output.trimmomatic_output}
    """

# ----------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------
rule fastqc_trimmed:
  input:
    trimmomatic_output= OUTPUTDIR + "02_trimmomatic/trimmomatic_output.txt"

  output:
    fastqc_trimmed_output= OUTPUTDIR + "03_fastqc/fastqc_trimmed_output.txt"

  params:
    fastqc_trimmed=expand( "05_Output/02_trimmomatic/{samples}.{run}.trimmed.fastq", samples=SAMPLES, run=RUN)

  conda: 
    CONTAINER + "fastqc.yaml"

  shell:
    """
    fastqc --outdir 05_Output/03_fastqc/ {params.fastqc_trimmed}
    touch {output.fastqc_trimmed_output}
    """

# ----------------------------------------------
# MultiQC to check the reads trimmed quality
# ----------------------------------------------

rule multiqc_trimmed:
  input:
    fastqc_trimmed_output= OUTPUTDIR + "03_fastqc/fastqc_trimmed_output.txt"

  output:
    multiqc_trimmed_output = OUTPUTDIR + "03_fastqc/trimmed_multiqc.txt"

  params:
    trim_qc = expand( "05_Output/03_fastqc/{samples}.{run}.trimmed_fastqc.zip", samples=SAMPLES, run=RUN),
    trim_multi_html = report(OUTPUTDIR + "03_fastqc/trimmed_multiqc.html", caption = ROOTDIR + "/report/multiqc.rst", category="01 quality report"), 
    multiqc_output_trim = OUTPUTDIR + "03_fastqc/trimmed_multiqc_data"

  conda:
    CONTAINER + "multiqc.yaml"

  shell: 
    """
    multiqc -n {params.trim_multi_html} {params.trim_qc} --force #run multiqc
    rm -rf {params.multiqc_output_trim} #clean-up
    touch {output.multiqc_trimmed_output}
    """
