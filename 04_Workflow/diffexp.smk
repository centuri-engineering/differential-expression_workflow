RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

# ----------------------------------------------
# DESeq2: Differential expression analysis
# ----------------------------------------------

def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
  output:
    rds = "05_Output/08_deseq2_init/all.rds"
    
  input:
    cts = "05_Output/07_cpm/count_filtered.txt",
    coldata = config["coldata"]
    
  conda:
    "../02_Container/deseq2.yaml"

  params:
    project = expand("{samples.project}",samples=samples.itertuples()),
    samples = expand("{samples.sample}",samples=samples.itertuples())

  threads: get_deseq2_threads()

  script:
    "../03_Script/deseq2.R"


rule diffexp:
  input:
    rds = "05_Output/08_deseq2_init/all.rds",
    rmd = "03_Script/diffexp.Rmd",
    pca = "03_Script/data_quality.R"
  
  output:
    #report("05_Output/09_pca/pca.pdf", caption="report/pca.rst", category="PCA")
    html_report = "05_Output/09_data_quality/diffexp.html"
    # Voir si figure en output
  
  params:
    pca_labels=config["pca"]["labels"]

  # singularity: 
  #   "docker://rocker/r-rmd"
  # conda:
  #   "../02_Container/deseq2.yaml"

  script:
    "../03_Script/diffexp_reports_compilation.R"