#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Workflow to benchmark and demonstrate Fastq compression with Squish
// Optionally compare results against Clumpify

// USAGE:
// nextflow run main.nf
// or
// nice nextflow run main.nf -profile docker --run_clumpify true --samplesheet samples.csv -resume

def splitList(value) {
    return (value ?: '')
        .toString()
        .split(/[;,]/)
        .collect { it.trim() }
        .findAll { it }
}


process RECOMPRESS_FASTQ {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_in), val(fastqout_base), path(paired_fastqs), val(paired_outputs)

    output:
    tuple val(sample_id), path(output_file), val(fastqout_base), path(paired_fastqs), val(paired_outputs)

    script:
    output_file = "recompress.${sample_id}.fastq.gz"
    """
    gunzip -c "${fastq_in}" | pigz -9 > "${output_file}"
    """
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
 *
 * When params.run_clumpify is true (default), Clumpify is also run on each
 * primary FASTQ and its results are included in report.csv for comparison.
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
      ${paired_arg} \\
      ${paired_out_arg} \\
      "${fastq_in}" \\
      "${fastqout_base}.${method}.fastq.gz"

    if [[ "${make_pdf}" == "true" ]]; then
      go tool pprof \\
        -pdf \\
        -output "\${method_outdir}/profile.${method}/memprofile.${method}.pdf" \\
        "\${method_outdir}/profile.${method}/mem.prof"
    fi

    test -s "\${method_outdir}/report.json"
    test -s "\${method_outdir}/manifest.txt"
    """
}

/*
 * Run Clumpify on a single primary FASTQ and emit a report.clumpify.json
 * whose fields match squish's report.json so both feed into the same CSV.
 *
 * uncompressed_bytes is null because Clumpify does not expose that count
 * without a separate decompression pass.
 */
process RUN_CLUMPIFY {
    tag "${sample_id}:clumpify"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_in), val(fastqout_base)
    val clump_k
    val clumpify_args

    output:
    path "${sample_id}/clumpify"                           , emit: clumpify_results
    path "${sample_id}/clumpify/report.clumpify.json"      , emit: report_json
    path "${sample_id}/clumpify/*.log"                     , emit: log

    script:
    def outdir    = "${sample_id}/clumpify"
    def out_fastq = "${fastqout_base}.clumpify.fastq.gz"
    def jvm_mem   = task.memory ? "-Xmx${(task.memory.toGiga() * 0.75).intValue()}g" : "-Xmx6g"
    """
    mkdir -p "${outdir}"

    input_size=\$(stat -L -c%s "${fastq_in}")

    clumpify.sh ${jvm_mem} \\
      in="${fastq_in}" \\
      out="${outdir}/${out_fastq}" \\
      k=${clump_k} \\
      tmpdir=./tmp \\
      ${clumpify_args} \\
      &> "${outdir}/${fastqout_base}.clumpify.log"

    output_size=\$(stat -c%s "${outdir}/${out_fastq}")
    size_diff=\$(( input_size - output_size ))
    size_ratio=\$(awk -v i="\${input_size}" -v o="\${output_size}" \\
      'BEGIN { if (i > 0) printf "%.6f", (i - o) / i; else print "0" }')

    cat > "${outdir}/report.clumpify.json" <<EOF
{
  "sort_method": "clumpify",
  "sort_engine": "clumpify",
  "clump_kmer_length": ${clump_k},
  "uncompressed_bytes": null,
  "output_size_bytes": \${output_size},
  "size_difference_bytes": \${size_diff},
  "size_reduction_ratio": \${size_ratio},
  "input": { "path": "${fastq_in}", "size_bytes": \${input_size} },
  "output": { "argument": "${out_fastq}" },
  "paired_outputs": []
}
EOF

    test -s "${outdir}/report.clumpify.json"
    """
}

/*
 * Convert a report JSON file (squish or Clumpify) into a single CSV row.
 * The sample_id is inferred from the grandparent directory name, which works
 * for both {sample_id}/{method}/report.json and
 * {sample_id}/clumpify/report.clumpify.json.
 */
def reportJsonToCsvRow(jsonFile) {
    def json        = new groovy.json.JsonSlurper().parseText(jsonFile.text)
    def sample_id   = jsonFile.parent.parent.name

    def header = [
        "sample_id",
        "sort_method",
        "sort_engine",
        "input_size_bytes",
        "uncompressed_bytes",
        "output_size_bytes",
        "size_difference_bytes",
        "size_reduction_ratio",
        "output_arg",
        "input_arg",
        "paired_output_args",
    ].join(",")

    def line = [
        sample_id,
        json["sort_method"],
        json["sort_engine"],
        json["input"]?.get("size_bytes") ?: "",
        json["uncompressed_bytes"] != null ? json["uncompressed_bytes"] : "",
        json["output_size_bytes"],
        json["size_difference_bytes"],
        json["size_reduction_ratio"],
        json["output"]["argument"],
        json["input"]["path"],
        (json["paired_outputs"] ?: []).collect { it["output"]["argument"] }.join(";"),
    ].join(",")

    return [header, line].join("\n")
}

workflow {
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def paired_fastqs = splitList(row.paired_fastqs).collect { file(it) }
            tuple(row.sample_id, file(row.fastqin), row.fastqout, paired_fastqs, row.paired_outputs ?: '')
        }

    if (params.recompress_r1 == true) {
        RECOMPRESS_FASTQ(samples_ch)
        samples_ch = RECOMPRESS_FASTQ.out
    }


    methods_ch = Channel
        .fromList(params.methods.tokenize(',').collect { it.trim() }.findAll { it })

    /*
     * combine creates the Cartesian product of samplesheet rows and sort
     * methods, so every sample is run through every configured sorter.
     */
    sample_methods_ch = samples_ch.combine(methods_ch)

    RUN_SQUISH_METHOD(
        sample_methods_ch,
        Channel.value(params.engine),
        Channel.value(params.bucket),
        Channel.value(params.buckets as int),
        Channel.value(params.clump_k as int),
        Channel.value(params.make_pdf as boolean)
    )

    all_reports_ch = RUN_SQUISH_METHOD.out.report_json

    if (params.run_clumpify) {
        clumpify_in_ch = samples_ch.map { sample_id, fastq_in, fastqout_base, _paired, _paired_out ->
            tuple(sample_id, fastq_in, fastqout_base)
        }

        RUN_CLUMPIFY(
            clumpify_in_ch,
            Channel.value(params.clump_k as int),
            Channel.value(params.clumpify_args)
        )

        all_reports_ch = all_reports_ch.mix(RUN_CLUMPIFY.out.report_json)
    }

    /*
     * Map every report JSON (squish and optionally Clumpify) to a CSV row
     * and collect into a single report.csv for side-by-side comparison.
     */
    all_reports_ch
        .map { jsonfile ->
            return reportJsonToCsvRow(jsonfile)
        }
        .collectFile(name: 'report.csv', storeDir: "${params.outdir}", newLine: true, keepHeader: true)
}
