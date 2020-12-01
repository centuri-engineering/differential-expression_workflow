RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
ref_level = config["diffexp"]["ref_level"]

# ----------------------------------------------
# DESeq2: Differential expression analysis
# ----------------------------------------------

def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
  input:
    cts = "05_Output/07_cpm/count_filtered.txt",
    #cts = "05_Output/07_cpm/cpm_filtered.txt",
    coldata = config["coldata"]

  output:
    rds = "05_Output/08_deseq2_init/all.rds"
        
  conda:
    "../02_Container/deseq2.yaml"

  params:
    project = expand("{samples.project}",samples=samples.itertuples()),
    samples = expand("{samples.sample}",samples=samples.itertuples()),
    ref_level = config["diffexp"]["ref_level"]

  threads: get_deseq2_threads()

  script:
    "../03_Script/deseq2.R"


rule diffexp:
  input:
    rds = "05_Output/08_deseq2_init/all.rds",
    rmd = "03_Script/diffexp.Rmd",
    pca = "03_Script/data_quality.R"
  
  output:
    html_report=report("05_Output/09_differential_expression/diffexp.html", caption="../report/diffexp.rst", category="03 Report differential expression"),
    table=report(expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}.csv", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="03 Report differential expression"),
    sur=report(expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}_signif-up-regulated.csv", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="03 Report differential expression"),
    sdr=report(expand("05_Output/09_differential_expression/{condition.condition}_vs_{ref_level}_signif-down-regulated.csv", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="03 Report differential expression")
  params:
    pca_labels=config["pca"]["labels"],
    ref_level = config["diffexp"]["ref_level"],
    lfcshrink_type = config["diffexp"]["lfcshrink_type"],
    gene_name = config["diffexp"]["gene_name"]

  script:
    "../03_Script/diffexp_reports_compilation.R"