##################################################################
# RNA velocity (kb nac + scVelo) and Scribe causal network
#
# Independent of SCEPTIC pseudotime: no rule here reads SCEPTIC output
# except compare_velocity_pseudotime, which only compares finished results.
# Velocity starts from results/integrated/integrated.h5ad (upstream of
# SCEPTIC), and Scribe genes are chosen from velocity fits, not pseudotime.
#
# Run (not added to rule all):
#   snakemake --executor slurm -j 32 velocity_all
##################################################################

NAC_DIR        = config.get("nac_ref_dir", "../reference/KB_Dmel_wMel_nac")
NAC_KALLISTO   = config.get("nac_kallisto", "/private/home/jomojaco/kallisto/build/src/kallisto")
VELOCITY_ENV   = config.get("velocity_env", "velocity")
SCRIBE_ENV     = config.get("scribe_env", "scribe")
SCRIBE_REPO    = config.get("scribe_repo", "../tools/Scribe-py")
VEL_SCRIPTS    = config.get("velocity_script_dir", "../snakemake_scripts/velocity")
VEL_PLATFORMS  = config.get("velocity_platforms", ["pipseq", "10x"])
SCRIBE_CHUNKS  = int(config.get("scribe_n_chunks", 10))
SCRIBE_PLATFORM = config.get("scribe_platform", "pipseq")
SCEPTIC_SAMPLE = config.get("velocity_sceptic_sample", "wolbachia_infection")

KB_TECH = {
    "pipseq": "0,0,8,0,11,17,0,20,26,0,31,39:0,39,51:1,0,0",
    "10x":    "10xv3",
}

def platform_of(sample_id):
    return str(samples_df.loc[sample_id][1]).strip()

def samples_for_platform(platform):
    return [s for s in SAMPLE_IDS if platform == "all" or platform_of(s) == platform]


rule velocity_all:
    input:
        expand("results/velocity/velocity_{platform}.h5ad", platform=VEL_PLATFORMS),
        "results/velocity/pseudotime_comparison/directionality_summary.csv",
        f"results/scribe/{SCRIBE_PLATFORM}/edges.csv",


# One-time: nascent + mature index from genome + GTF (Dmel + wMel combined)
rule build_nac_index:
    input:
        genome = config.get("nac_genome_fasta", "MISSING_nac_genome_fasta"),
        gtf    = config.get("nac_gtf", "MISSING_nac_gtf"),
    output:
        idx     = f"{NAC_DIR}/index.idx",
        t2g     = f"{NAC_DIR}/t2g.txt",
        cdna    = f"{NAC_DIR}/cdna.fa",
        nascent = f"{NAC_DIR}/nascent.fa",
        c1      = f"{NAC_DIR}/cdna.t2c",
        c2      = f"{NAC_DIR}/nascent.t2c",
    log: "logs/velocity/build_nac_index.log"
    threads: 16
    resources:
        slurm_partition = "medium",
        mem_mb  = 128000,
        runtime = hms_to_minutes("4:00:00"),
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate kallisto_bustools; set -u
        kb ref --workflow=nac --kallisto {NAC_KALLISTO} -t {threads} \
            -i {output.idx} -g {output.t2g} \
            -c1 {output.c1} -c2 {output.c2} \
            -f1 {output.cdna} -f2 {output.nascent} \
            {input.genome} {input.gtf}
        """


# Parallel to map_pipseq / map_10x; existing count matrices are untouched
rule map_nac:
    input:
        read1 = lambda wc: get_fastq_files(wc.sample_id)[0],
        read2 = lambda wc: get_fastq_files(wc.sample_id)[1],
        idx   = rules.build_nac_index.output.idx,
        t2g   = rules.build_nac_index.output.t2g,
        c1    = rules.build_nac_index.output.c1,
        c2    = rules.build_nac_index.output.c2,
    output:
        h5ad = "results/nac/{sample_id}/adata_nac.h5ad",
    params:
        outdir = "results/nac/{sample_id}",
        tech   = lambda wc: KB_TECH[platform_of(wc.sample_id)],
    log: "logs/velocity/map_nac_{sample_id}.log"
    threads: config["pipseeker_threads"]
    resources:
        slurm_partition = config["pipseeker_partition"],
        mem_mb  = config.get("nac_mem", 128000),
        runtime = hms_to_minutes(config.get("nac_time", "6:00:00")),
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate kallisto_bustools; set -u
        # Remove partial kb output from a failed run; kb count fails if the folder exists
        rm -rf {params.outdir}
        kb count --workflow=nac --kallisto {NAC_KALLISTO} \
            -i {input.idx} -g {input.t2g} -c1 {input.c1} -c2 {input.c2} \
            -x {params.tech} -o {params.outdir} -t {threads} --h5ad \
            {input.read1} {input.read2}
        mv {params.outdir}/counts_unfiltered/adata.h5ad {output.h5ad}
        """


rule rna_velocity:
    input:
        integrated = "results/integrated/integrated.h5ad",
        nac = lambda wc: expand("results/nac/{sample_id}/adata_nac.h5ad",
                                sample_id=samples_for_platform(wc.platform)),
    output:
        h5ad = "results/velocity/velocity_{platform}.h5ad",
        qc   = "results/velocity/{platform}/unspliced_fraction_by_sample.csv",
        fits = "results/velocity/{platform}/velocity_gene_fits.csv",
    wildcard_constraints:
        platform = "pipseq|10x|all"
    params:
        fig_dir      = "results/velocity/{platform}",
        neighbor_rep = config.get("velocity_neighbor_rep", "X_pca_harmony"),
    log: "logs/velocity/rna_velocity_{platform}.log"
    threads: 16
    resources:
        slurm_partition = "medium",
        mem_mb  = config.get("velocity_mem", 256000),
        runtime = hms_to_minutes(config.get("velocity_time", "12:00:00")),
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate {VELOCITY_ENV}; set -u
        python {VEL_SCRIPTS}/run_velocity.py \
            --integrated {input.integrated} \
            --nac-files {input.nac} \
            --platform {wildcards.platform} \
            --neighbor-rep {params.neighbor_rep} \
            --out-h5ad {output.h5ad} \
            --fig-dir {params.fig_dir} \
            --n-jobs {threads}
        """


# The only rule that sees both results
rule compare_velocity_pseudotime:
    input:
        velocity = "results/velocity/velocity_pipseq.h5ad",
        sceptic  = f"results/sceptic/{SCEPTIC_SAMPLE}/sceptic_{SCEPTIC_SAMPLE}.h5ad",
    output:
        summary = "results/velocity/pseudotime_comparison/directionality_summary.csv",
        rank    = "results/velocity/pseudotime_comparison/rank_agreement.csv",
    params:
        outdir = "results/velocity/pseudotime_comparison",
    log: "logs/velocity/compare_velocity_pseudotime.log"
    threads: 4
    resources:
        slurm_partition = "medium",
        mem_mb  = 64000,
        runtime = hms_to_minutes("2:00:00"),
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate {VELOCITY_ENV}; set -u
        python {VEL_SCRIPTS}/compare_velocity_pseudotime.py \
            --velocity {input.velocity} --sceptic {input.sceptic} --outdir {params.outdir}
        """


rule scribe_select_genes:
    input:
        velocity   = "results/velocity/velocity_{platform}.h5ad",
        regulators = config.get("scribe_regulator_list", "MISSING_scribe_regulator_list"),
    output:
        regs    = "results/scribe/{platform}/regulators.txt",
        targets = "results/scribe/{platform}/targets.txt",
    params:
        outdir    = "results/scribe/{platform}",
        n_targets = config.get("scribe_n_targets", 300),
        max_regs  = config.get("scribe_max_regulators", 100),
        exclude   = (f"--exclude {config['scribe_exclude_genes']}"
                     if config.get("scribe_exclude_genes") else ""),
    log: "logs/scribe/select_{platform}.log"
    threads: 1
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = 30,
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate {SCRIBE_ENV}; set -u
        export PYTHONPATH={SCRIBE_REPO}:${{PYTHONPATH:-}}
        python {VEL_SCRIPTS}/run_scribe.py --mode select \
            --velocity {input.velocity} --regulators {input.regulators} {params.exclude} \
            --n-targets {params.n_targets} --max-regulators {params.max_regs} \
            --outdir {params.outdir}
        """


rule scribe_chunk:
    input:
        velocity = "results/velocity/velocity_{platform}.h5ad",
        regs     = "results/scribe/{platform}/regulators.txt",
    output:
        "results/scribe/{platform}/chunks/rdi_chunk{chunk}.csv",
    wildcard_constraints:
        chunk = r"\d{3}"
    params:
        outdir    = "results/scribe/{platform}",
        max_cells = config.get("scribe_max_cells", 5000),
    log: "logs/scribe/{platform}_chunk{chunk}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb  = 32000,
        runtime = hms_to_minutes(config.get("scribe_chunk_time", "12:00:00")),
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate {SCRIBE_ENV}; set -u
        export PYTHONPATH={SCRIBE_REPO}:${{PYTHONPATH:-}}
        python {VEL_SCRIPTS}/run_scribe.py --mode chunk \
            --velocity {input.velocity} --outdir {params.outdir} \
            --chunk $((10#{wildcards.chunk})) --n-chunks {SCRIBE_CHUNKS} \
            --max-cells {params.max_cells}
        """


rule scribe_merge:
    input:
        expand("results/scribe/{{platform}}/chunks/rdi_chunk{chunk:03d}.csv",
               chunk=range(SCRIBE_CHUNKS)),
    output:
        edges = "results/scribe/{platform}/edges.csv",
        rdi   = "results/scribe/{platform}/rdi_matrix.csv",
    params:
        outdir = "results/scribe/{platform}",
    log: "logs/scribe/{platform}_merge.log"
    threads: 1
    resources:
        slurm_partition = "medium", mem_mb = 16000, runtime = 30,
    shell:
        """
        exec > {log} 2>&1
        set +u; source {CONDA_BASE}/etc/profile.d/conda.sh; conda activate {SCRIBE_ENV}; set -u
        python {VEL_SCRIPTS}/run_scribe.py --mode merge --outdir {params.outdir}
        """
