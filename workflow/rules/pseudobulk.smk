rule paired_pseudobulk:
    input:
        audit=rules.audit_inputs.output.manifest,
        script="workflow/scripts/pseudobulk.R"
    output:
        counts="results/pseudobulk/counts.tsv",
        groups="results/pseudobulk/groups.tsv",
        design="results/pseudobulk/design.tsv",
        tested="results/pseudobulk/feature_filter.tsv",
        da="results/pseudobulk/differential_accessibility.tsv",
        plots="results/pseudobulk/diagnostics.pdf",
        summary="results/pseudobulk/summary.json",
        session="results/pseudobulk/sessionInfo.txt"
    conda: "../envs/pseudobulk.yaml"
    log: "logs/pseudobulk/model.log"
    benchmark: "benchmarks/pseudobulk/model.tsv"
    resources: mem_mb=16000
    shell:
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla {input.script:q} "
        "{input.audit:q} > {log:q} 2>&1"


rule pseudobulk_report:
    input:
        summary=rules.paired_pseudobulk.output.summary,
        provenance=rules.audit_inputs.output.manifest,
        script="workflow/scripts/report.py"
    output: TARGET
    log: "logs/pseudobulk/report.log"
    resources: mem_mb=256
    shell:
        "{PYTHON:q} {input.script:q} --title 'Paired donor-level differential accessibility' "
        "--summary {input.summary:q} --provenance {input.provenance:q} --output {output:q} > {log:q} 2>&1"
