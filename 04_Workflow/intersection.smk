# ----------------------------------------------
# Presence or absence of differential expressed genes in different comparisons
# ----------------------------------------------

rule intersection:
  input:
    diffexp_output = OUTPUTDIR + "09_differential_expression/diffexp.txt"

  output:
    intersection_output = OUTPUTDIR + "10_intersection/intersection.txt"

  params:
    gene_list_1 = OUTPUTDIR + "09_differential_expression/" + config["intersection"]["gene_list_1"],
    gene_list_2 = OUTPUTDIR + "09_differential_expression/" + config["intersection"]["gene_list_2"],
    common_genes = report(OUTPUTDIR + "10_intersection/" + config["intersection"]["gene_list_1"] + "_vs_" + config["intersection"]["gene_list_2"] + ".txt", category="03 Report differential expression"),
    stat_file_1 = OUTPUTDIR + "09_differential_expression/" + config["intersection"]["stat_file_1"],
    stat_file_2 = OUTPUTDIR + "09_differential_expression/" + config["intersection"]["stat_file_2"],
    stat_common_genes = report(OUTPUTDIR + "10_intersection/common_genes_stat_" + config["intersection"]["gene_list_1"] + "_vs_" + config["intersection"]["gene_list_2"] + ".txt", category="03 Report differential expression"),

  conda:
    CONTAINER + "intersection.yaml"

  script:
    SCRIPTDIR + "intersection.R"
