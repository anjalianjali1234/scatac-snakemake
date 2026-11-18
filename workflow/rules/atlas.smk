# All biological jobs use the same versioned R environment.
rule atlas_qc:
    input:
        audit=rules.audit_inputs.output.manifest,
        script="workflow/scripts/atlas.R"
    output:
        obj="results/atlas/qc.rds",
        metrics="results/atlas/qc_metrics.csv",
        figure="results/atlas/qc.pdf",
        session="results/atlas/session_qc.txt"
    conda: "../envs/atac.yaml"
    log: "logs/atlas/qc.log"
    benchmark: "benchmarks/atlas/qc.tsv"
    resources: mem_mb=16000
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla {input.script:q} qc "
        "{input.audit:q} > {log:q} 2>&1"


rule atlas_lsi:
    input:
        obj=rules.atlas_qc.output.obj,
        audit=rules.audit_inputs.output.manifest,
        script="workflow/scripts/atlas.R"
    output:
        obj="results/atlas/lsi.rds",
        decisions="results/atlas/cell_filtering.csv",
        correlations="results/atlas/depth_correlations.csv",
        figure="results/atlas/lsi_depth.pdf",
        session="results/atlas/session_lsi.txt"
    conda: "../envs/atac.yaml"
    log: "logs/atlas/lsi.log"
    benchmark: "benchmarks/atlas/lsi.tsv"
    resources: mem_mb=16000
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla {input.script:q} lsi "
        "{input.audit:q} > {log:q} 2>&1"


rule atlas_embedding:
    input:
        obj=rules.atlas_lsi.output.obj,
        audit=rules.audit_inputs.output.manifest,
        script="workflow/scripts/atlas.R"
    output:
        obj="results/atlas/atlas.rds",
        cells="results/atlas/cells.csv",
        figures="results/atlas/embedding_and_activity.pdf",
        summary="results/atlas/summary.json",
        session="results/atlas/session_embedding.txt"
    conda: "../envs/atac.yaml"
    log: "logs/atlas/embedding.log"
    benchmark: "benchmarks/atlas/embedding.tsv"
    resources: mem_mb=16000
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla {input.script:q} embedding "
        "{input.audit:q} > {log:q} 2>&1"


rule atlas_report:
    input:
        summary=rules.atlas_embedding.output.summary,
        provenance=rules.audit_inputs.output.manifest,
        script="workflow/scripts/report.py"
    output: TARGET
    log: "logs/atlas/report.log"
    resources: mem_mb=256
    shell:
        "{PYTHON:q} {input.script:q} --title 'PBMC accessibility atlas: exploratory analysis' "
        "--summary {input.summary:q} --provenance {input.provenance:q} --output {output:q} > {log:q} 2>&1"
