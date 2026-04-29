#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

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
    path "${method}/report*json", emit : report_json

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

    RUN_SQUISH_METHOD.out.report_json.map{ jsonFile ->
        def jsonSlurper = new JsonSlurper()
        def Map json = (Map) new JsonSlurper().parseText(jsonFile.text)
        def sort_method = json["sort_method"]
        def sort_engine = json["sort_engine"]
        def uncompressed_bytes = json["uncompressed_bytes"]
        def output_size_bytes = json["output_size_bytes"]
        def size_difference_bytes = json["size_difference_bytes"]
        def size_reduction_ratio = json["size_reduction_ratio"]

        return [
            sort_method,
            sort_engine,
            uncompressed_bytes,
            output_size_bytes,
            size_difference_bytes,
            size_reduction_ratio
            ].join(",")
    }.collectFile(name: 'report.csv', storeDir: "${params.outdir}", newLine: true)
}
