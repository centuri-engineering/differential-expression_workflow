SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
# rule fastqc:
#   output:
#     expand( "05_Output/01_fastqc/{samples}_{ext}_fastqc.html", samples=SAMPLES, ext=["1","2"]),
#     expand( "05_Output/01_fastqc/{samples}_{ext}_fastqc.zip", samples=SAMPLES, ext=["1","2"])

#   input:
#     expand( "00_RawData/{samples}_{ext}.fastq", samples=SAMPLES, ext=["1","2"])

#   conda: 
#     "../02_Container/fastqc.yaml"

#   shell:
#     "fastqc --outdir 05_Output/01_fastqc/ {input}"

# ----------------------------------------------
# Trimmomatic: trimming reads and removing adapter sequences
# ----------------------------------------------
rule trimmomatic:
  input:
    sample=expand("00_RawData/{samples}_{ext}.fastq", samples=SAMPLES, ext=["1","2"])

  output:
    sample_trimmed=expand( "05_Output/02_trimmomatic/{samples}_{ext}.trimmed.fastq", samples=SAMPLES, ext=["1","2"]),
    sample_untrimmed=expand( "05_Output/02_trimmomatic/{samples}_{ext}un.trimmed.fastq", samples=SAMPLES, ext=["1","2"])

  conda: 
    "../02_Container/trimmomatic.yaml"

  shell:
    """
    sample=({input.sample})
    sample_trimmed=({output.sample_trimmed})
    sample_untrimmed=({output.sample_untrimmed})
    len=${{#sample[@]}}
    for (( i=0; i<$len; i=i+2 ))
    do trimmomatic PE -threads 4 ${{sample[$i]}} ${{sample[$i+1]}} ${{sample_trimmed[$i]}} ${{sample_untrimmed[$i]}} ${{sample_trimmed[$i+1]}} ${{sample_untrimmed[$i+1]}} LEADING:20 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
    done
    """

# ----------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------
rule fastqc_trimmed:
  output:
    report(expand( "05_Output/03_fastqc/{samples}_{ext}.trimmed_fastqc.html", samples=SAMPLES, ext=["1","2"]), caption="../report/fastqc.rst", category="01 Row data quality"),
    expand( "05_Output/03_fastqc/{samples}_{ext}.trimmed_fastqc.zip", samples=SAMPLES, ext=["1","2"])

  input:
    expand( "05_Output/02_trimmomatic/{samples}_{ext}.trimmed.fastq", samples=SAMPLES, ext=["1","2"])

  conda: 
    "../02_Container/fastqc.yaml"

  shell:
    "fastqc --outdir 05_Output/03_fastqc/ {input}"