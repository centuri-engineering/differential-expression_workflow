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
    cpm_filtering_output = OUTPUTDIR + "07_cpm/cpm_filtering.txt"

  output:
    deseq2_init_output = OUTPUTDIR + "08_deseq2_init/deseq2_init.txt"

  conda:
    "../02_Container/deseq2.yaml"

  params:
    cts = OUTPUTDIR + "07_cpm/count_filtered.txt",
    coldata = config["coldata"],
    rds = OUTPUTDIR + "08_deseq2_init/all.rds",
    normalized_counts_file = report(OUTPUTDIR + "08_deseq2_init/normalized_counts.tsv", caption="../report/normalized_counts.rst", category="02 Count matrices"),
    project = expand("{samples.project}",samples=samples.itertuples()),
    samples = expand("{samples.sample}",samples=samples.itertuples()),
    ref_level = config["diffexp"]["ref_level"],
    rmproj_list = config["filtering"]["rmproj"]

  threads: get_deseq2_threads()

  script:
    "../03_Script/deseq2.R"


rule diffexp:
  input:
    deseq2_init_output = OUTPUTDIR + "08_deseq2_init/deseq2_init.txt"
  output:
    diffexp_output = OUTPUTDIR + "09_differential_expression/diffexp.txt"
  conda:
    CONTAINER + "diffexp.yaml"
  params:
    rds = "05_Output/08_deseq2_init/all.rds",
    rmd = "03_Script/diffexp.Rmd",
    pca = "03_Script/data_quality.R",
    coldata = config["coldata"],
    html_report=report(OUTPUTDIR + "09_differential_expression/diffexp.html", caption="../report/diffexp.rst", category="03 Report differential expression"),
    pca_labels=config["pca"]["labels"],
    ref_level = config["diffexp"]["ref_level"],
    lfcshrink_type = config["diffexp"]["lfcshrink_type"],
    gene_name = config["diffexp"]["gene_name"],
    mutant_level = config["diffexp"]["mutant_level"],
    nbpval = config["diffexp"]["nbpval"],
    FCcutoff=config["diffexp"]["FCcutoff"],
    pCutoff=config["diffexp"]["pCutoff"],
  message: 
    "Run the differential expression analyses"  
  script:
    SCRIPTDIR + "diffexp_reports_compilation.R"
