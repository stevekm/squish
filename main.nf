#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Nextflow version of the Makefile test-run-all recipes.
 *
 * The workflow runs the same five methods used by the Makefile:
 *
 *   alpha, gc, qual, alpha-heap, clump
 *
 * Use --engine memory to mirror `make test-run-all`.
 * Use --engine external to mirror `make test-run-all-external`.
 *
 * The squish binary is expected to be available on PATH. For local runs, place
 * it in a local bin directory that Nextflow prepends to PATH.
 */

process RUN_SQUISH_METHOD {
    tag "${method}:${engine}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    val method
    path fastq_in
    val engine
    val bucket
    val buckets
    val fastqout_base
    val make_pdf

    output:
    path "${method}", emit: method_results

    /*
     * Each method gets its own result directory. This avoids filename
     * collisions between memory and external runs and keeps the profile,
     * order, PDF, and output FASTQ together.
     */
    script:
    """
    method_outdir="${method}"
    mkdir -p "\${method_outdir}"

    echo ">>> RUNNING: method=${method} engine=${engine} bucket=${bucket} buckets=${buckets}"

    squish \\
      -outdir "\${method_outdir}" \\
      -engine "${engine}" \\
      -bucket "${bucket}" \\
      -buckets "${buckets}" \\
      -m "${method}" \\
      -orderFile "order.${method}.txt" \\
      -reportFile "report.${method}.json" \\
      -memProf "mem.${method}.prof" \\
      -cpuProf "cpu.${method}.prof" \\
      "${fastq_in}" \\
      "${fastqout_base}.${method}.fastq.gz"

    if [[ "${make_pdf}" == "true" ]]; then
      go tool pprof \\
        -pdf \\
        -output "\${method_outdir}/profile.${method}/memprofile.${method}.pdf" \\
        "\${method_outdir}/profile.${method}/mem.${method}.prof"
    fi

    test -s "\${method_outdir}/report.${method}.json"
    """
}

workflow {
    fastq_in = file(params.fastqin)

    methods_ch = Channel
        .fromList(params.methods.tokenize(',').collect { it.trim() }.findAll { it })

    RUN_SQUISH_METHOD(
        methods_ch,
        Channel.value(fastq_in),
        Channel.value(params.engine),
        Channel.value(params.bucket),
        Channel.value(params.buckets as int),
        Channel.value(params.fastqout),
        Channel.value(params.make_pdf as boolean)
    )
}
