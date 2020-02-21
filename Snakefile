localrules: generate_min_p, gen_gwsig_snps, get_block_sizes, pull_genotypes_condsig_position, all_pull_genotypes_condsig_position, merge_cond_results, cond_final_compare, condsig_create_vcf, annot_condsig, run_final_cond, annot_gwsig, create_gwsig_vcf, create_igv_gwsig
import os

def read_file_array(fname):
    with open(fname) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 
    return(content)

TRAITS = read_file_array(os.environ['PHENO_NAMES'])

# to begin with the sex chromosomes exist separately as XY and X - we later merge them into chr 24
CHROMS_somatic = list(range(1,23))
CHROMS = ["XY", "X"]
CHROMS += CHROMS_somatic

# use this once the sex chromosomes have been collapsed into CHR 24
CHROMS_num = CHROMS_somatic + ['23']
## function which returns a list of blocks for a chromsome?
def get_blocks(Wildcards):
    return(True)

rule process_boltlmm:
    input:
    output:
        "data/all_variants/{trait}.assoc",
        "data/assoc_files/{trait}_gwsig.assoc"
    params:
        section="base"
    resources:
        mem_mb=lambda wildcards, attempt: ((attempt * 10000)+80000)
    benchmark: "benchmarks/process_boltlmm_all_variants_{trait}.tsv"
    shell: "Rscript ./process_boltlmm_all_variants.R {wildcards.trait} "

rule initial_mapping_file_from_bolt:
    input: expand("data/assoc_files/{trait}_gwsig.assoc", trait=TRAITS)
    output: "data/tables/mapping_file.vcf"
    params: section="base"
    resources: mem_mb=1000
    shell: "Rscript ./inital_mapping_file_from_bolt.R"


rule snp_stats_rsids:
    input: "data/tables/mapping_file.vcf"
    output:
        "data/tables/mapping_file.vcf.gz",
        "data/tables/mapping_file_rsids.vcf"
    params: section="base"
    resources: mem_mb=3000
    shell: "./snp_stats_rsids.sh"

rule create_mapping_rsid:
    input: "data/tables/mapping_file_rsids.vcf"
    output:
        "data/tables/mapping_file_rsid.csv",
        "data/tables/mapping_file_all.csv"
    params: section="base"
    resources: mem_mb=1000
    shell: "Rscript ./create_mapping_rsid.R"

rule all_process_boltlmm:
    input: expand("data/assoc_files/{trait}_gwsig.assoc", trait=TRAITS)
    output: "data/all_process_boltlmm.txt"
    params:
        section="base"
    resources:
        mem_mb=100
    shell: "echo all done > data/all_process_boltlmm.txt"

rule generate_min_p:
    input: expand("data/assoc_files/{trait}_gwsig.assoc", trait=TRAITS)
    output: "data/all_min.tsv"
    params:
        section="base"
    resources:
        mem_mb=11000
    shell: "Rscript ./generate_min_p.R"

rule gen_gwsig_snps:
    input: "data/all_min.tsv"
    output: expand("data/genfiles/gwsig_snps_chr{chromosome}_pos.tsv", chromosome=CHROMS)
    params:
        section="base"
    resources:
        mem_mb=2000
    shell: "Rscript ./gen_gwsig_snps.R"

rule block_traits:
    input: expand("data/{trait}_gwsig.assoc", trait=TRAITS)
    output: dynamic("data/blocks/pull_ids_block_{blockid}.tsv")
    params:
        section="base"
    resources:
        mem_mb=12000
    threads: 10
    shell: "Rscript ./block_traits.R"

rule transform_block_list:
    input: dynamic("data/blocks/pull_ids_block_{blockid}.tsv")
    output:
        "data/blocks/block_chrs_order.tsv",
        "data/blocks/blocks_each_chrom.tsv"
    params: section="base"
    resources: mem_mb=500
    shell: "Rscript ./transform_blocks_list.R"

rule get_block_sizes:
    input:
        "data/blocks/block_chrs_order.tsv",
        "data/blocks/blocks_each_chrom.tsv",
        "data/blocks/pull_ids_block_1.tsv"
    output: "data/blocks/block_sizes_db.csv"
    params: section="base"
    resources: mem_mb=700
    shell: """
    source activate itr
    python ./get_block_sizes.py
    """

rule pull_genotypes_gwsig:
    input: "data/genfiles/gwsig_snps_chr{chromosome}_pos.tsv"
    output:
        "data/genfiles/chr{chromosome}_gwsig_clean_id.gen",
        "data/genfiles/ids_{chromosome}.tsv"
    params:
        section="base"
    resources:
        mem_mb=lambda wildcards, attempt: ((attempt * 5000)+5000)
    threads: 2
    shell: "pull_genotypes_gwsig.sh {wildcards.chromosome}"

rule sex_chrom_create_new_samples:
    input:
    output:
        "data/sex_chrom_x.sample",
        "data/sex_chrom_xy.sample"
    params: section="base"
    resources: mem_mb=12000
    shell: "Rscript ./sex_chrom_create_new_samples.R"

rule sex_chrom_merge_gen_files:
    input:
        "data/genfiles/chrX_gwsig_clean_id.gen", # comes from pull_genotypes_gwsig
        "data/genfiles/chrXY_gwsig_clean_id.gen",
        "data/sex_chrom_x.sample",
        "data/sex_chrom_xy.sample"
    output:
        "data/genfiles/chrXY_merged.sample",
        "data/genfiles/chrXY_merged.gen"
    params: section="base"
    resources: mem_mb=5000
    shell: "./sex_chrom_merge_gen_files.sh"

rule all_pull_genotypes_gwsig:
    input: expand("data/genfiles/ids_{chromosome}.tsv", chromosome=CHROMS)
    output: "data/all_pull_genotypes_gwsig.txt"
    params:
        section="base"
    resources:
        mem_mb=100
    shell: "echo aaaa > data/all_pull_genotypes_gwsig.txt"


# this rule is fake and shouldn't be run, you should run the script manually
rule job_creation_script:
    input:
        expand("data/genfiles/ids_{chromosome}.tsv", chromosome=CHROMS_somatic),
        expand("data/genfiles/chr{chromosome}_gwsig_clean_id.gen", chromosome=CHROMS_somatic),
        "data/genfiles/chrXY_merged.gen",
        "data/genfiles/chrXY_merged.sample",
        "data/blocks/pull_ids_block_1.tsv"
    output: "data/condout/blocks/ind_var_block1.tsv"
    params: section="post_block"
    resources: mem_mb=200
    shell: "python ./job_creation_script.py"

rule merge_univar_results:
    input: "data/condout/blocks/ind_var_block1.tsv"
    output: expand("data/condout/trait/var_list_{trait}.tsv", trait=TRAITS)
    params: section="post_block"
    resources: mem_mb=500
    shell: "Rscript ./merge_univar_results.R"
    
rule table_genotypes_condsig:
    input: expand("data/condout/trait/var_list_{trait}.tsv", trait=TRAITS)
    output: expand("data/condout/pull_ids/pull_ids_chr{chrom}_pos.tsv", chrom=CHROMS_num)
    params: section="post_block"
    resources: mem_mb=500
    shell: "Rscript ./table_genotypes_condsig.R"

rule pull_genotypes_condsig_position:
    input: "data/condout/pull_ids/pull_ids_chr{chrom}_pos.tsv"
    output:
        "data/condout/subgenhd5/ids_pos_{chrom}.tsv",
        "data/condout/subgenhd5/chr{chrom}_pos_clean_id.gen",
        "data/condout/subgenhd5/chr{chrom}_pos_clean_id.h5"
    params: section="post_block"
    resources: mem_mb=6000
    shell: "./pull_genotypes_condsig_position.sh {wildcards.chrom}"

rule all_pull_genotypes_condsig_position:
    input:
        expand("data/condout/subgenhd5/chr{chrom}_pos_clean_id.gen", chrom=CHROMS_num),
        expand("data/condout/subgenhd5/chr{chrom}_pos_clean_id.h5", chrom=CHROMS_num)
    output: "data/all_genotype_pulls.txt"
    params: section="post_block"
    resources: mem_mb=100
    shell: "echo aaaaa > data/all_genotype_pulls.txt"

rule run_final_cond:
    input:
        "data/condout/subgenhd5/chr{chrom}_pos_clean_id.gen",
        "data/condout/subgenhd5/chr{chrom}_pos_clean_id.h5"
    output: "data/condout/results_chrom/condsig_{trait}_gwas_normalised_chr_{chrom}.tsv"
    params: section="post_block"
    resources: mem_mb=250000
    shell: "Rscript ./cond.R {wildcards.trait} {wildcards.chrom}"
        

rule merge_cond_results:
    input: expand("data/condout/results_chrom/condsig_{{trait}}_gwas_normalised_chr_{chrom}.tsv", chrom=CHROMS_num)
    output: "data/condout/results/condsig_{trait}_gwas_normalised.tsv"
    params: section="post_block"
    resources: mem_mb=100
    shell: "Rscript ./merge_cond_results.R {wildcards.trait}"

rule cond_final_compare:
    input: expand("data/condout/results/condsig_{trait}_gwas_normalised.tsv", trait=TRAITS)
    output: "data/ukbb_compare_previous_number_signals.csv"
    params: section="post_block"
    resources: mem_mb=200
    shell: "Rscript ./cond_final_compare.R"

rule condsig_create_vcf:
    input: expand("data/condout/results/condsig_{trait}_gwas_normalised.tsv", trait=TRAITS)
    output: "data/condout/vep/all_condsig.vcf"
    params: section="post_block"
    resources: mem_mb=1000
    shell: "Rscript ./condsig_create_vcf.R"

rule annot_condsig:
    input: "data/condout/vep/all_condsig.vcf"
    output: "data/condout/vep/tabular_all_condsig.txt"
    params: section="post_block"
    resources: mem_mb=1000
    shell: "./annot_condsig.sh"

rule create_gwsig_vcf:
    input: "data/assoc_files/{trait}_gwsig.assoc"
    output: "data/vcf/{trait}_gwsig.vcf"
    params: section="not_cond"
    resources: mem_mb=3000
    shell: "Rscript ./create_gwsig_vcf.R {wildcards.trait}"

rule annot_gwsig:
    input: "data/vcf/{trait}_gwsig.vcf"
    output: "data/vep/{trait}_gwsig_vep.vcf"
    params: section="not_cond"
    resources: mem_mb=3000
    shell: "./annot.sh {wildcards.trait}"

rule create_igv_gwsig:
    input:
        expand("data/vep/{trait}_gwsig_vep.vcf", trait=TRAITS),
        "data/tables/mapping_file_rsid.csv"
    output: expand("data/condout/igv/{trait}_gwsig.gwas", trait=TRAITS)
    params: section="not_cond"
    resources: mem_mb=3000
    shell: "Rscript ./create_igv.R gwsig"

#        ancient("data/condout/subgenhd5/ids_pos_{chrom}.tsv"),
#        ancient("data/condout/subgenhd5/chr{chrom}_pos_clean_id.h5")
rule calculate_ld_dosage:
    input:
    output: "data/condout/ldclump_dosage/dosage_{chrom}.ld"
    params: section="not_cond"
    resources: mem_mb=3000
    conda: "python_envs/itr.yaml"
    shell: "python ./calculate_ld_dosage.py --chrom {wildcards.chrom}"

rule merge_calculate_ld_dosage:
    input: expand("data/condout/ldclump_dosage/dosage_{chrom}.ld", chrom=CHROMS_num)
    output: "data/condout/ldclump_dosage/dosage.ld"
    params: section="not_cond"
    resources: mem_mb=300
    shell: "python ./merge_calculate_ld_dosage.py"

rule clumping_algorithm:
    input: "data/condout/ldclump_dosage/dosage.ld"
    output: "data/condout/ldclump_dosage/hg19_clump_table.tsv"
    params: section="not_cond"
    resources: mem_mb=1000
    shell: "Rscript ./clumping_algorithm.R"

rule ld_clump:
    input: expand("data/condout/results/condsig_{trait}_gwas_normalised.tsv", trait=TRAITS)
    output: "data/condout/ldclump/hg19_clump_table.tsv"
    params: section="not_cond"
    resources: mem_mb=3000
    shell: "./ld_clump.sh"

rule create_igv_cond:
    input:
        expand("data/condout/results/condsig_{trait}_gwas_normalised.tsv", trait=TRAITS),
        expand("data/condout/assoc_files/{trait}_gwsig.assoc", trait=TRAITS),
        "data/condout/ldclump_dosage/hg19_clump_table.tsv",
        "data/condout/vep/tabular_all_condsig.txt"
    output:
        "data/condout/igv_orig_names/{trait}_condsig.gwas"
    params: section="not_cond"
    resources: mem_mb=200
    shell: "Rscript ./create_igv.R cond"

rule create_final_results:
    input:
        "data/condout/vep/tabular_all_condsig.txt",
        "data/condout/ldclump/hg19_clump_table.tsv",
        expand("data/condout/igv_orig_names/{trait}_condsig.gwas", trait=TRAITS)
    output:
        "data/BCX_final_output_raw.csv",
        "data/BCX_final_output.tsv"
    resources: mem_mb=2000
    params: section="not_cond"
    shell: "Rscript ./create_final_results.R"

rule astle_compare_pull_table:
    input:
    output: expand("data/astle_compare/pull_tables/{chrom}.csv", chrom=CHROMS_somatic)
    resources: mem_mb=100
    params: section="not_cond"
    shell: "Rscript ./astle_compare_pull_table.R"

rule astle_compare_pull_genotypes:
    input: "data/astle_compare/pull_tables/{chrom}.csv"
    output:
        "data/astle_compare/genfiles/ids_{chrom}.tsv",
        "data/astle_compare/genfiles/chr{chrom}_clean_id.gen"
    resources: mem_mb=500
    params: section="not_cond"
    shell: "./astle_compare_pull_genotypes.sh {wildcards.chrom}"

rule astle_compare_calculate_ld:
    input: "data/astle_compare/genfiles/chr{chrom}_clean_id.gen"
    output: "data/astle_compare/ld/ld_chr{chrom}.ld"
    resources: mem_mb=2000
    conda: "python_envs/itr.yaml"
    params: section="not_cond"
    shell: "python ./astle_compare_calculate_ld.py --chrom {wildcards.chrom}"

rule astle_compare_find_novel_clumps:
    input:
        expand("data/astle_compare/ld/ld_chr{chrom}.ld", chrom=CHROMS_somatic),
        "data/condout/ldclump_dosage/h19_clump_table.tsv"
    output: "data/astle_compare/hg19_clump_table_novel.tsv"
    resources: mem_mb=5000
    conda: "python_envs/itr.yaml"
    params: section="not_cond"
    shell: "python ./astle_compare_find_novel_clumps.py"

rule locuszoom_input:
    input: "data/all_variants/{trait}.assoc"
    output: "data/locuszoom/input/{trait}.tsv"
    resources: mem_mb=5000
    params: section="not_cond"
    shell: "Rscript ./locuszoom_input.R {wildcards.trait}"

rule locuszoom_input_all:
    input: expand("data/locuszoom/input/{trait}.tsv", trait=TRAITS)
    output: "data/locuszoom/input/all_done.txt"
    resources: mem_mb=200
    params: section="not_cond"
    shell: "echo aaa > data/locuszoom/input/all_done.txt"

rule locuszoom_wrapper:
    input: expand("data/locuszoom/input/{trait}.tsv", trait=TRAITS)
    output: "data/locuszoom/output/done.txt"
    resources: mem_mb=2000
    params: section="not_cond"
    shell: "./locuszoom_script_wrapper.sh"
