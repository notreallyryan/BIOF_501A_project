#!/usr/bin/env nextflow
/*

    A basic pipeline for core genome phylogenetic analysis of paired end illumina V. cholerae samples.

    biof_501a_project:
        1. download samples fron NCBI SRA (if necessary)
        2. download BUSCO and CARD databases (if necessary)
        3. QC and trim raw reads using fastp
        4. Short read assembly using SKESA
        5. assembly QC using BUSCO and Quast
        6. core genome phylogeny using PARSNP
        7. Antimicrobial resistance gene identification using CARD RGI
        8. Plot tree and AMR heatmap using ggtree in R.
        9. Generate MULTIQC report.

    Usage:
        nextflow run main.nf -profile docker \
            --input_dir /path/to/input_dir \
            --input_csv /path/to/samplesheet.csv \
            --outdir results

    CSV columns:
        accession  (required)  SRA run accession, e.g. SRR0000001
        fastq1     (optional)  path relative to --input_dir
        fastq2     (optional)  path relative to --input_dir; leave blank for single-end
        country    (required)  country of origin
        year       (required)  year of sampling

    Samples whose FASTQ file(s) are missing from --input_dir are downloaded from NCBI SRA.
*/

nextflow.enable.dsl = 2

/*
    --------------------------------------------
    MODULES / SUBWORKFLOWS
    --------------------------------------------
*/
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS } from './subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main'
include { FASTP   } from './modules/nf-core/fastp/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SKESA } from './modules/nf-core/skesa/main'
include { BUSCO_DOWNLOAD } from './modules/nf-core/busco/download/main'
include { BUSCO_BUSCO } from './modules/nf-core/busco/busco/main'
include { QUAST } from './modules/nf-core/quast/main'
include { RGI_DOWNLOADDB } from './modules/ebi-metagenomics/rgi/downloaddb/main'
include { RGI_CARDANNOTATION } from './modules/nf-core/rgi/cardannotation/main'
include { RGI_MAIN } from './modules/nf-core/rgi/main/main'
include { PARSNP } from './modules/nf-core/parsnp/main'
include { R_PLOT } from './modules/local/r_plots/main'


/*
    --------------------------------------------
    MAIN WORKFLOW
    --------------------------------------------
*/
workflow {

    // Validate inputs 
    if (!params.input_dir) { error "Missing required parameter: --input_dir" }
    if (!params.input_csv) { error "Missing required parameter: --input_csv" }

    def input_dir = file(params.input_dir)
    if (!input_dir.isDirectory()) {
        error "Input directory does not exist or is not a directory: ${params.input_dir}"
    }

    /*
       --------------------------------------------
       1.  DOWNLOAD SAMPLE DATA IF NECESSARY
       --------------------------------------------
    */

    // Read the CSV and split into local vs. to-download
    ch_samples = Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.containsKey('accession')) {
                error "Input CSV is missing required column: accession"
            }
            def accession = row.accession?.trim()
            if (!accession) {
                error "Input CSV has a row with an empty accession: ${row}"
            }

            def fq1_rel = row.fastq1?.trim()
            if (!accession) {
                error "Input CSV has a row with an empty fastq1: ${row}"
            }

            def fq2_rel = row.fastq2?.trim()
            if (!accession) {
                error "Input CSV has a row with an empty fastq2: ${row}"
            }

            def fq1     = fq1_rel ? input_dir.resolve(fq1_rel) : null
            def fq2     = fq2_rel ? input_dir.resolve(fq2_rel) : null

            def present = (fq1 && fq1.exists() && fq2 && fq2.exists())
            if (!present) {
                log.warn "FASTQ files for ${accession} not found in ${params.input_dir}; will download from NCBI SRA."
            }

            def reads = [ fq1, fq2 ].findAll { it } 
            [ accession, reads, present ]
        }
        .branch { accession, reads, present ->
            local:    present
                return [ accession, reads ]
            download: true
                return accession
        }

    // Download the missing samples 
    // Subworkflow input: [ val(meta), val(sra_id) ], plus an optional dbGaP key (not needed for this)
    ch_sra_ids = ch_samples.download.map { accession -> [ [ id: accession ], accession ] }
    FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS( ch_sra_ids, [] )

    // Downloaded reads: [ meta, reads ] -> [ accession, [ reads ] ]
    ch_downloaded = FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
        .map { meta, reads ->
            def fastqs = (reads instanceof List ? reads : [ reads ]).sort { it.name }
            [ meta.id, fastqs ]
        }

    // Combine into one channel and add the meta map the modules expect
    ch_reads = ch_samples.local
        .mix(ch_downloaded)
        .map { accession, reads ->
            [ [ id: accession, single_end: reads.size() == 1 ], reads ]
        }


    /* 
       --------------------------------------------
       2.  DOWNLOAD DATABASES IF NECESSARY
       --------------------------------------------
    */

    // BUSCO DB
    def buscoDownloads = file("${params.databases}/busco_downloads")
    def buscoLineage   = file("${buscoDownloads}/lineages/${params.BUSCO_db}")

    // check if already exists. If not, then download.
    def ch_busco_db
    if (buscoLineage.isDirectory()) {
        log.info "BUSCO database found in ${buscoDownloads}."
        ch_busco_db = Channel.value(buscoDownloads)
    } else {
        log.info "BUSCO database not found in ${params.databases}. Downloading now."
        // need the specific species for best results, stored in nextflow.config under BUSCO_db
        BUSCO_DOWNLOAD(params.BUSCO_db)
        ch_busco_db = BUSCO_DOWNLOAD.out.download_dir
    }

    // CARD DB
    def cardDownloads = file("${params.databases}/card_dir")
    def card_json = file("${cardDownloads}/card.json")

    // check if already exists. If not, then download.
    def ch_card_db
    if (card_json.exists()) {
        log.info "CARD database found in ${cardDownloads}."
        ch_card_db = Channel.value(cardDownloads)
    } else {
        log.info "CARD database not found in ${params.databases}. Downloading now."
        RGI_DOWNLOADDB()
        ch_card_db = RGI_DOWNLOADDB.out.card_json
    }


    /*
        --------------------------------------------
        3. FASTP FOR RAW READ QC AND TRIM
        --------------------------------------------
    */
    
    // FASTP input: 
    // [ val(meta), path(reads), path(adapter_fasta) ], (adapter_fasta is not necessary for fastp)
    // bool(discard_trimmed_pass),
    // bool(save_trimmed_fail),
    // bool(save_merged)

    ch_fastp_input = ch_reads .map { meta, reads -> [meta, reads, []]}
    FASTP(ch_fastp_input, false, false, false)

    /*
       --------------------------------------------
       4.  ASSEMBLY
       --------------------------------------------
    */

    // SKESA input:
    // [ val(meta), path(reads) ]
    ch_skesa_input = FASTP.out.reads
    SKESA(ch_skesa_input)

    /*
       --------------------------------------------
       5.  ASSEMBLY QC
       --------------------------------------------
    */

    // BUSCO input:
    // [ val(meta), path(reads) ]
    // val(mode) (set to genome mode)
    // val lineage (the species to download if database is not downloaded yet.)
    // path busco_lineages_path 
    // path config_file
    // val clean_intermediates

    BUSCO_BUSCO(
        SKESA.out.fasta,
        'genome',
        params.BUSCO_db,
        ch_busco_db,
        [],
        true
    )

    // QUAST input:
    // tuple val(meta) , path(consensus) - assemblies to evaluate
    // tuple val(meta2), path(fasta) - reference genome (optional)
    // tuple val(meta3), path(gff) - reference annotation (optional)

    // pull up the v cholera reference genome and gff file
    def reference   = file("${params.ref_genome}")
    def gff   = file("${params.ref_gff}")

    ch_ref_input = [ [ id: 'V_cho_reference' ], reference ]
    ch_gff_input = [ [ id: 'V_cho_gff' ], gff ]

    QUAST(SKESA.out.fasta, ch_ref_input, ch_gff_input)

    /*
       --------------------------------------------
       6.  PHYLOGENETIC ANALYSIS
       --------------------------------------------
    */

    ch_panaroo_input = SKESA.out.fasta
        .map { meta, fasta -> fasta }
        .collect()
        .map { fasta -> [ [ id: 'parsnp' ], fasta ] }

    // PARSNP input:
    // [ val(meta), path(assemblies) ]
    // path(reference)
    PARSNP(ch_panaroo_input, reference)


    /*
       --------------------------------------------
       7.  AMR PROFILING
       --------------------------------------------
    */
    // Preprocess CARD database
    RGI_CARDANNOTATION(ch_card_db)

    // running CARD
    RGI_MAIN(
        SKESA.out.fasta,
        ch_card_db,
        []
    )

    /*
       --------------------------------------------
       8.  R PLOTTING
       --------------------------------------------
    */
    ch_metadata = Channel.fromPath(params.input_csv) 
    ch_tree = PARSNP.out.tree.map { meta, tree -> tree }
    ch_rgi = RGI_MAIN.out.tsv.map { meta, tsv -> tsv }.collect()

    R_PLOT(ch_tree, ch_metadata, ch_rgi)

    /*
       --------------------------------------------
       9.  MULTIQC REPORT
       --------------------------------------------
    */

    ch_multiqc_input = Channel.empty()
        .mix(FASTP.out.json.map { meta, json -> json })
        .mix(BUSCO_BUSCO.out.short_summaries_txt.map { meta, summary -> summary })
        .mix(QUAST.out.results.map {meta, results -> results} )
        .collect()
        .map { files -> [ [ id: 'multiqc' ], files, [], [], [], [] ] }

    MULTIQC(ch_multiqc_input)
}
