#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

def splitList(value) {
    return (value ?: '')
        .toString()
        .split(/[;,]/)
        .collect { it.trim() }
        .findAll { it }
}

/*
 * Nextflow version of the Makefile test-run-all recipes.
 *
 * The workflow runs the same four methods used by the Makefile:
 *
 *   alpha, gc, qual, clump
 *
 * Use --engine memory to mirror `make test-run-all`.
 * Use --engine external to mirror `make test-run-all-external`.
 *
 * Input FASTQ files are read from a CSV samplesheet with these columns:
 *
 *   sample_id,fastqin,fastqout,paired_fastqs,paired_outputs
 *
 * paired_fastqs and paired_outputs are optional semicolon-separated lists.
 * paired_outputs are filename bases; the workflow appends .<method>.fastq.gz.
 *
 * The squish binary is expected to be available on PATH. For local runs, place
 * it in a local bin directory that Nextflow prepends to PATH.
 */

process RUN_SQUISH_METHOD {
    tag "${sample_id}:${method}:${engine}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_in), val(fastqout_base), path(paired_fastqs), val(paired_outputs), val(method)
    val engine
    val bucket
    val buckets
    val clump_k
    val make_pdf

    output:
    path "${sample_id}/${method}", emit: method_results
    path "${sample_id}/${method}/report*json", emit : report_json
    path "${sample_id}/${method}/manifest*txt", emit : manifest_txt

    /*
     * Each sample/method pair gets its own result directory. This avoids
     * filename collisions and keeps the profile, order, PDF, report, and
     * output FASTQ together.
     */
    script:
    def paired_files = paired_fastqs instanceof List ? paired_fastqs : (paired_fastqs ? [paired_fastqs] : [])
    def paired_arg = paired_files ? "-paired \"${paired_files.join(',')}\"" : ""
    def paired_output_items = splitList(paired_outputs).collect { "${it}.${method}.fastq.gz" }
    def paired_out_arg = paired_output_items ? "-pairedOut \"${paired_output_items.join(',')}\"" : ""
    """
    method_outdir="${sample_id}/${method}"
    mkdir -p "\${method_outdir}"

    echo ">>> RUNNING: sample=${sample_id} method=${method} engine=${engine} bucket=${bucket} buckets=${buckets} clumpK=${clump_k}"

    squish \\
      -outdir "\${method_outdir}" \\
      -engine "${engine}" \\
      -bucket "${bucket}" \\
      -buckets "${buckets}" \\
      -clumpK "${clump_k}" \\
      -m "${method}" \\
      -orderFile "order.${method}.txt" \\
      -reportFile "report.${method}.json" \\
      -manifestFile "manifest.${method}.txt" \\
      -memProf "mem.${method}.prof" \\
      -cpuProf "cpu.${method}.prof" \\
      ${paired_arg} \\
      ${paired_out_arg} \\
      "${fastq_in}" \\
      "${fastqout_base}.${method}.fastq.gz"

    if [[ "${make_pdf}" == "true" ]]; then
      go tool pprof \\
        -pdf \\
        -output "\${method_outdir}/profile.${method}/memprofile.${method}.pdf" \\
        "\${method_outdir}/profile.${method}/mem.${method}.prof"
    fi

    test -s "\${method_outdir}/report.${method}.json"
    test -s "\${method_outdir}/manifest.${method}.txt"
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def paired_fastqs = splitList(row.paired_fastqs).collect { file(it) }
            tuple(row.sample_id, file(row.fastqin), row.fastqout, paired_fastqs, row.paired_outputs ?: '')
        }

    methods_ch = Channel
        .fromList(params.methods.tokenize(',').collect { it.trim() }.findAll { it })

    /*
     * combine creates the Cartesian product of samplesheet rows and sort
     * methods, so every sample is run through every configured sorter.
     */
    sample_methods_ch = samples_ch
        .combine(methods_ch)

    RUN_SQUISH_METHOD(
        sample_methods_ch,
        Channel.value(params.engine),
        Channel.value(params.bucket),
        Channel.value(params.buckets as int),
        Channel.value(params.clump_k as int),
        Channel.value(params.make_pdf as boolean)
    )

    /*
     * Each process emits one JSON report. Convert those reports into compact
     * CSV rows for quick comparison across samples and methods.
     */
    RUN_SQUISH_METHOD.out.report_json.map{ jsonFile ->
        def jsonSlurper = new JsonSlurper()
        def Map json = (Map) new JsonSlurper().parseText(jsonFile.text)
        def sort_method = json["sort_method"]
        def sort_engine = json["sort_engine"]
        def sample_id = jsonFile.parent.parent.name
        def uncompressed_bytes = json["uncompressed_bytes"]
        def output_size_bytes = json["output_size_bytes"]
        def size_difference_bytes = json["size_difference_bytes"]
        def size_reduction_ratio = json["size_reduction_ratio"]
        def output_arg = json["output"]["argument"]
        def input_arg = json["input"]["path"]
        def paired_output_args = (json["paired_outputs"] ?: []).collect { it["output"]["argument"] }.join(";")

        header = [
            "sample_id",
            "sort_method",
            "sort_engine",
            "uncompressed_bytes",
            "output_size_bytes",
            "size_difference_bytes",
            "size_reduction_ratio",
            "output_arg",
            "input_arg",
            "paired_output_args"
            ].join(",")
        line = [
            sample_id,
            sort_method,
            sort_engine,
            uncompressed_bytes,
            output_size_bytes,
            size_difference_bytes,
            size_reduction_ratio,
            output_arg,
            input_arg,
            paired_output_args
            ].join(",")
        output = [header, line].join("\n")

        return output
    }.collectFile(name: 'report.csv', storeDir: "${params.outdir}", newLine: true, keepHeader: true)
}
